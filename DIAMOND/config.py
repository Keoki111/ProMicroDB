# ================= 配置文件 =================
# 注意：请优先修改 config.py 中的参数，常规使用无需修改 main.py 的核心逻辑。

import os # 这个不要动！

# [路径设置]
# 建议：INPUT_FOLDER 尽量只放待处理的fasta文件、config.py和main.py，避免后缀干扰。
INPUT_FOLDER = "./input_fasta"                          # 输入文件夹 (不是具体文件)。 
OUTPUT_FILE  = "./output/DME_result.tsv" # 结果保存路径。
# 结果格式是tsv，如果写成了csv则excel打开会乱码）。

# [文件类型支持]
FILE_EXTENSIONS = ["*.fa", "*.fasta", "*.faa", "*.pep", "*.txt"]

# [Diamond 运行参数]
# fast: 默认 (极快) | sensitive: 灵敏 | very-sensitive: 非常灵敏 | ultra-sensitive: 很慢但很准
# 请合理设置线程数！
DIAMOND_MODE = "ultra-sensitive"
THREADS      = 100        # 请先查看现有进程数和cpu占用，避免设置太高影响他人，实在不确定建议使用 nice 命令
E_VALUE      = "1e-10"    # E-value
MIN_IDENTITY = 30.0      # 最小相似度 (0-100)
MIN_COVERAGE = 80.0       # 最小覆盖度 (0-100)，如果不设覆盖度限制则填入“ 0.0 ”

# [输出数量控制]
# 设置每个查询序列 (Query) 保留多少条匹配结果。
# 本质上你输入0和一个很大的数字（比如说10000）没有区别。
# 原则是可多不可少，比如你输入了100，实际只有30个符合要求的结果，那最终会输出30个，不会报错。
# 0  = 输出【所有】符合阈值 (E-value/Identity/Coverage) 的结果 (不做数量截断)。
# 1  = 只输出最佳匹配 (Best Hit)。
# n  = n为任意正整数，代表输出前 n 个匹配。
MAX_HITS = 0

# [输出内容模式]
# 1 = 【原生模式】仅输出 Diamond 原始比对数据，格式为Output Format 6 (qseqid, sseqid, pident, evalue...)。
#     -> 特点：非常快，不需要加载物种数据库，内存占用不大。
# 2 = 【标准模式】部分原生数据 + 匹配的物种分类信息 (合并后)。
#     -> 特点：信息最全，包含相似度、覆盖度以及对应的物种分类信息。
# 3 = 【精简模式】仅输出 ID (Query, Subject, Genome) + 物种分类信息。
#     -> 特点：文件体积小，去除了比对分数、长度等冗余列，只留下物种分类信息。
OUTPUT_STYLE = 2

# [高级功能：物种分类处理]
# 仅在 OUTPUT_STYLE 为 2 或 3 时生效。
# True  = 拆分 Lineage 为 Domain...Species 多列，并清洗前缀。
# False = 保持 GTDB 长字符串格式。
SPLIT_TAXONOMY = True

# 这里强烈建议参考数据库和元数据表都放在一个独立且固定的位置。
# 这样只需要设定一次，以后都不需要再修改了，非常省事。
DIAMOND_DB   = "./database/uhgp-95.dmnd"
TAX_FILE     = "./database/genomes-all_metadata.tsv"