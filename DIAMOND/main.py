# 这是主程序

import pandas as pd
import subprocess
import os
import glob
import sys
from io import StringIO
import config

# ==================== 核心功能函数 ====================

def clean_taxonomy(df):
    # 这个函数用于清洗并拆分 GTDB 物种层级信息

    if not config.SPLIT_TAXONOMY:
        return df

    print("  -> 正在执行物种拆分与清洗...", flush=True)
    
    target_col = 'Lineage'
    # 自动寻找 Lineage 列
    if target_col not in df.columns:
        if len(df.columns) > 1:
            target_col = df.columns[1] 
        else:
            return df

    # 1. 拆分 Lineage
    ranks = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
    # 使用 expand=True 直接拆分成多列
    split_df = df[target_col].str.split(';', expand=True)
    
    # 防止层级不足或过多
    if split_df.shape[1] > len(ranks):
        split_df = split_df.iloc[:, :len(ranks)]
    split_df.columns = ranks[:split_df.shape[1]]
    
    # 2. 清洗前缀 (d__, p__ 等)
    prefixes = {
        'Domain': '^d__', 'Phylum': '^p__', 'Class': '^c__', 
        'Order': '^o__',  'Family': '^f__', 'Genus': '^g__', 'Species': '^s__'
    }
    for col in split_df.columns:
        if col in prefixes:
            split_df[col] = split_df[col].str.replace(prefixes[col], '', regex=True).str.strip()

    # 3. 合并回原表 (去掉旧的 Lineage 列)
    df_new = pd.concat([df.drop(columns=[target_col]), split_df], axis=1)
    return df_new

def load_taxonomy_map():
    # 这个函数用于加载物种注释元数据表
    # 注意：只有当 OUTPUT_STYLE 为 2 或 3 时才会被调用

    print(f"正在加载物种分类数据库: {config.TAX_FILE} ...", flush=True)
    try:
        df = pd.read_csv(config.TAX_FILE, sep='\t')
        
        # 规范化 Genome_ID
        df.rename(columns={df.columns[0]: 'Genome_ID'}, inplace=True)
        df['Genome_ID'] = df['Genome_ID'].astype(str).str.strip()
        
        if config.SPLIT_TAXONOMY:
            df = clean_taxonomy(df)
            
        return df
    except Exception as e:
        print(f"Error: 加载物种文件失败: {e}")
        sys.exit(1)

def run_pipeline(fasta_path, df_tax=None):
    # 对单个文件运行 Diamond 并根据配置处理结果

    # 1. 设置 Diamond 的 hits 数量
    # Diamond 中 --max-target-seqs 0 代表无限制输出 (report all alignments)
    # 也可以设为非常大的数字来模拟 'all'
    k_hits = "0" if config.MAX_HITS == 0 else str(config.MAX_HITS)
    
    cmd = [
        "diamond", "blastp",
        "--db", config.DIAMOND_DB,
        "--query", fasta_path,
        "--outfmt", "6", "qseqid", "sseqid", "pident", "length", "evalue", "qlen", "slen",
        "--threads", str(config.THREADS),
        "--evalue", config.E_VALUE,
        "--max-target-seqs", k_hits, 
        "--quiet"
    ]
    
    if config.DIAMOND_MODE != "fast":
        cmd.append(f"--{config.DIAMOND_MODE}")

    # 2. 执行并捕获输出
    result = subprocess.run(cmd, capture_output=True, text=True)
    if not result.stdout:
        return None 

    # 3. 解析为 DataFrame
    cols = ["qseqid", "sseqid", "pident", "length", "evalue", "qlen", "slen"]
    df = pd.read_csv(StringIO(result.stdout), sep='\t', names=cols)
    
    # 提取 Genome_ID (用于后续关联)
    df['sseqid'] = df['sseqid'].astype(str).str.strip()
    df['Genome_ID'] = df['sseqid'].str.rsplit('_',n=1).str[0]
    
    # 计算 Coverage
    max_len = df[['qlen', 'slen']].max(axis=1)
    df['Coverage'] = (df['length'] / max_len) * 100
    
    # 4. 阈值过滤
    mask_ident = df['pident'] >= config.MIN_IDENTITY
    mask_cov = df['Coverage'] >= config.MIN_COVERAGE
    df_clean = df[mask_ident & mask_cov].copy()
    
    if df_clean.empty:
        return None

    # 5. 添加样本名列
    df_clean.insert(0, 'Sample', os.path.basename(fasta_path))
    
    # ================= 分支逻辑：根据 OUTPUT_STYLE 格式化输出 =================
    
    # 【模式 1】：原生模式 - 直接返回，不合并物种
    if config.OUTPUT_STYLE == 1:
        # 此时 df_tax 是 None，也不需要操作
        return df_clean

    # 此时如果 df_tax 为 None (代码逻辑错误防护)，则报错
    if df_tax is None:
        raise ValueError("需要物种数据库但未加载 (df_tax is None)")

    # 【合并物种信息】 (Left Join)
    # 使用 Genome_ID 进行连接
    df_merged = df_clean.merge(df_tax, on='Genome_ID', how='left')

    # 【模式 3】：精简模式 - 剔除比对统计列
    if config.OUTPUT_STYLE == 3:
        # 定义要移除的比对统计列
        stats_cols_to_drop = ['pident', 'length', 'evalue', 'qlen', 'slen', 'Coverage']
        # 注意：保留 Sample, qseqid, sseqid, Genome_ID 以及所有来自 df_tax 的列
        df_final = df_merged.drop(columns=stats_cols_to_drop, errors='ignore')
        return df_final
    
    # 【模式 2】：标准模式 - 返回全集 (默认)
    # 删除中间计算的长度列，保持整洁，但保留 pident, evalue, coverage 等
    df_final = df_merged.drop(columns=['qlen', 'slen'], errors='ignore')
    
    return df_final

# ==================== 主程序 ====================

def main():
    print(f"="*60)
    print(f"UHGP Pipeline Started")
    print(f"Mode: {config.DIAMOND_MODE} | Threads: {config.THREADS}")
    print(f"Hits per Query: {'ALL' if config.MAX_HITS == 0 else config.MAX_HITS}")
    print(f"Output Style: {config.OUTPUT_STYLE} (1=Raw, 2=Std, 3=Slim)")
    print(f"="*60)
    
    if not os.path.exists(config.DIAMOND_DB):
        print(f"Error: 数据库不存在: {config.DIAMOND_DB}")
        sys.exit(1)
        
    # 优化：只有当模式不是 1 时，才加载巨大的物种表
    df_tax = None
    if config.OUTPUT_STYLE in [2, 3]:
        df_tax = load_taxonomy_map()
    else:
        print("模式为 1 (Raw)，跳过物种数据库加载，直接输出比对结果。")
    
    # 扫描文件
    input_files = []
    for ext in config.FILE_EXTENSIONS:
        pattern = os.path.join(config.INPUT_FOLDER, ext)
        input_files.extend(glob.glob(pattern))
    input_files = sorted(list(set(input_files)))

    if not input_files:
        print(f"在 {config.INPUT_FOLDER} 没找到文件。")
        sys.exit(1)

    print(f"检测到 {len(input_files)} 个样本。")
    
    all_results = []
    
    for i, f in enumerate(input_files):
        print(f"[{i+1}/{len(input_files)}] Processing: {os.path.basename(f)} ... ", end="", flush=True)
        try:
            res = run_pipeline(f, df_tax)
            if res is not None:
                all_results.append(res)
                print(f"Hits: {len(res)}")
            else:
                print("No Hits")
        except Exception as e:
            print(f"\n[Warning] 处理文件 {os.path.basename(f)} 时出错: {e}")
            # 单个文件出错不中断整个流程
            continue
            
    if all_results:
        print(f"正在合并所有结果并保存...")
        final_table = pd.concat(all_results, ignore_index=True)
        
        # 保存
        final_table.to_csv(config.OUTPUT_FILE, sep='\t', index=False)
        print(f"\n成功！")
        print(f"结果已保存至: {config.OUTPUT_FILE}")
        print(f"总行数: {len(final_table)}")
    else:
        print("\n流程结束，未产生任何有效结果。")

if __name__ == "__main__":
    main()