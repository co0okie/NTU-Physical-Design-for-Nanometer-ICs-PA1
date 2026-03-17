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

    # 2. 依序生成 Nets
    for i in range(1, num_nets + 1):
        # 隨機決定這條 net 要連幾個 cell (2, 3 或 4 個，各 33% 機率)
        num_connected_cells = random.choice([2, 3, 4])
        
        # 防呆機制：如果總 cell 數量比要抽的數量還少
        num_connected_cells = min(num_connected_cells, num_cells)

        # random.sample 保證抽出來的 cell 不會重複
        selected_cells = random.sample(cells_pool, num_connected_cells)

        # 組合字串並印出 (加上 PA1 規定的結尾分號)
        cells_str = " ".join(selected_cells)
        print(f"NET n{i} {cells_str} ;")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="PA1 F-M Partitioning 測資產生器")
    parser.add_argument("-n", "--nets", type=int, required=True, help="Net 的數量 (必填)")
    parser.add_argument("-c", "--cells", type=int, required=True, help="Cell 的數量 (必填)")
    parser.add_argument("-b", "--balance", type=float, help="Balance degree (選填，0~1之間)")

    args = parser.parse_args()
    generate_testcase(args.nets, args.cells, args.balance)