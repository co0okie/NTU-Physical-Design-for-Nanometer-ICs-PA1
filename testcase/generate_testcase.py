#! /usr/bin/python3

import argparse
import random

def generate_testcase(num_nets, num_cells, balance_degree=None):
    # 1. 決定 Balance Degree
    if balance_degree is None:
        # 使用常態分佈，平均值設為 0.2，標準差設為 0.15
        # 透過 while 迴圈確保數值絕對落在 0 ~ 1 之間
        while True:
            b = random.gauss(0.2, 0.15)
            if 0 < b < 1.0:
                balance_degree = round(b, 2)
                break
    else:
        if not (0 <= balance_degree <= 1):
            raise ValueError("Balance degree 必須介於 0 到 1 之間")

    # 印出 Balance Degree
    print(f"{balance_degree}")

    # 建立 Cell 的候選池 (c1, c2, ..., cn)
    cells_pool = [f"c{i}" for i in range(1, num_cells + 1)]

    # 設定長尾分佈的形狀參數 (alpha)
    # alpha 越小 (例如 1.1)，尾巴越長，越容易出現極端巨大的 Net
    # alpha 越大 (例如 2.5)，數值會越集中在底線附近 (2, 3, 4...)
    pareto_alpha = 1.5 

    # 2. 依序生成 Nets
    for i in range(1, num_nets + 1):
        # random.paretovariate(alpha) 會產生 >= 1.0 的浮點數
        # 我們取整數後加 1，這樣最小值必定是 2 (當 raw_val 落在 1.0 ~ 1.99 之間時)
        raw_val = random.paretovariate(pareto_alpha)
        num_connected_cells = int(raw_val) + 1
        
        # 防呆機制：確保單一 Net 連接的 Cell 數量不會超過全場總 Cell 數
        num_connected_cells = min(num_connected_cells, num_cells)

        # random.sample 保證抽出來的 cell 不會重複
        selected_cells = random.sample(cells_pool, num_connected_cells)

        # 組合字串並印出 (加上 PA1 規定的結尾分號)
        cells_str = " ".join(selected_cells)
        print(f"NET n{i} {cells_str} ;")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PA1 F-M Partitioning 測資產生器 (長尾分佈版)")
    parser.add_argument("-n", "--nets", type=int, required=True, help="Net 的數量 (必填)")
    parser.add_argument("-c", "--cells", type=int, required=True, help="Cell 的數量 (必填)")
    parser.add_argument("-b", "--balance", type=float, help="Balance degree (選填，0~1之間)")

    args = parser.parse_args()
    generate_testcase(args.nets, args.cells, args.balance)