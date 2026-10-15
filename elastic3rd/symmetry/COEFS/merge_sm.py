import os
from monty.serialization import dumpfn
def read_strains_to_voigt(filepath):
    """
    从文件中读取多个应变矩阵（3x3），返回 nx6 的 Voigt 向量列表。
    使用逐行读取 + 空行检测，兼容各种换行符和多余空白。
    """
    voigt_list = []
    current_matrix = []

    with open(filepath, 'r') as f:
        for line in f:
            stripped = line.strip()
            # 如果是空行（包括只有空白字符的行），且当前正在收集矩阵
            if not stripped:
                if current_matrix:          # 遇到空行，结束当前矩阵
                    if len(current_matrix) != 3:
                        raise ValueError(f"矩阵行数应为3，实际得到{len(current_matrix)}行")
                    # 将 current_matrix 转为 Voigt 向量
                    e11 = current_matrix[0][0]
                    e22 = current_matrix[1][1]
                    e33 = current_matrix[2][2]
                    e23 = 2 * current_matrix[1][2]
                    e13 = 2 * current_matrix[0][2]
                    e12 = 2 * current_matrix[0][1]
                    voigt_list.append([e11, e22, e33, e23, e13, e12])
                    current_matrix = []      # 重置
                # 如果 current_matrix 为空，则跳过连续空行
                continue

            # 非空行：解析数字
            parts = stripped.split()
            if len(parts) != 3:
                raise ValueError(f"每行应有3个数字，实际得到: {stripped}")
            row = list(map(float, parts))
            current_matrix.append(row)

    # 处理文件末尾可能没有空行的情况
    if current_matrix:
        if len(current_matrix) != 3:
            raise ValueError(f"矩阵行数应为3，实际得到{len(current_matrix)}行")
        e11 = current_matrix[0][0]
        e22 = current_matrix[1][1]
        e33 = current_matrix[2][2]
        e23 = 2 * current_matrix[1][2]
        e13 = 2 * current_matrix[0][2]
        e12 = 2 * current_matrix[0][1]
        voigt_list.append([e11, e22, e33, e23, e13, e12])

    return voigt_list

# 使用示例
if __name__ == "__main__":
    files = ['C1', 'C2', 'H1', 'H2', 'T1', 'T2', 'R1', 'R2', 'O', 'M', 'N']
    order = 3
    for file in files:
        results = {}
        for mode in ['ENERGY', 'STRESS2']:
            filename = os.path.join(str(order), mode, file)
            result_lst = read_strains_to_voigt(filename)
            results[mode.capitalize()] = result_lst
        dumpfn(results, f'{file}.json', indent=4)
