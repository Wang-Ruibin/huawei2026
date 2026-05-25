import os
import numpy as np
import matplotlib.pyplot as plt
import random


# ===================== 读取输入 =====================
def read_input(file_path):
    with open(file_path, 'r') as f:
        lines = f.read().strip().splitlines()

    poly1 = np.array(list(map(float, lines[1].split()))).reshape(-1, 2)
    poly2 = np.array(list(map(float, lines[2].split()))).reshape(-1, 2)

    moves = []
    for line in lines[4:]:
        x, y = map(float, line.split())
        moves.append((x, y))

    return poly1, poly2, moves


# ===================== 读取 MSV =====================
def read_output(file_path):
    msv = []
    with open(file_path, 'r') as f:
        for line in f:
            x, y = map(float, line.split())
            msv.append((x, y))
    return msv


# ===================== 初始图 =====================
def plot_initial(poly1, poly2, save_dir):
    plt.figure(figsize=(6, 6))

    p1 = np.vstack([poly1, poly1[0]])
    p2 = np.vstack([poly2, poly2[0]])

    plt.plot(p1[:, 0], p1[:, 1], 'b-', label='Polygon 1')
    plt.plot(p2[:, 0], p2[:, 1], 'r-', label='Polygon 2')

    plt.title("Initial Position")
    plt.legend()
    plt.axis('equal')
    plt.grid(True)

    plt.savefig(os.path.join(save_dir, "initial.png"))
    plt.close()


# ===================== 对比图 =====================
def plot_compare(poly1, poly2, move, sep_vec, step, move_idx, save_dir):
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))

    # ===== 左：移动后 =====
    ax = axes[0]

    p1 = np.vstack([poly1, poly1[0]])
    moved_poly2 = poly2 + np.array(move)
    p2 = np.vstack([moved_poly2, moved_poly2[0]])

    ax.plot(p1[:, 0], p1[:, 1], 'b-')
    ax.plot(p2[:, 0], p2[:, 1], 'r-')

    centroid = poly2.mean(axis=0)
    ax.arrow(centroid[0], centroid[1],
             move[0], move[1],
             head_width=0.02,
             length_includes_head=True,
             color='orange')

    ax.set_title(f"Moved (#{move_idx})")
    ax.axis('equal')
    ax.grid(True)

    # ===== 右：分离后 =====
    ax = axes[1]

    ax.plot(p1[:, 0], p1[:, 1], 'b-')
    ax.plot(p2[:, 0], p2[:, 1], 'r--')  # 原位置（虚线）

    separated_poly2 = moved_poly2 + np.array(sep_vec)
    p3 = np.vstack([separated_poly2, separated_poly2[0]])
    ax.plot(p3[:, 0], p3[:, 1], 'g-')

    moved_centroid = centroid + np.array(move)
    ax.arrow(moved_centroid[0], moved_centroid[1],
             sep_vec[0], sep_vec[1],
             head_width=0.02,
             length_includes_head=True,
             color='purple')

    ax.set_title("After Separation")
    ax.axis('equal')
    ax.grid(True)

    # ===== 保存（简洁命名）=====
    filename = os.path.join(save_dir, f"s{step}_m{move_idx}.png")
    plt.savefig(filename)
    plt.close()


# ===================== 主流程 =====================
random.seed(42)  # 可复现

for i in range(1, 8):
    input_file = f"practice_{i}.in"
    output_file = f"practice_{i}.out"
    save_dir = f"practice_{i}"

    os.makedirs(save_dir, exist_ok=True)

    poly1, poly2, moves = read_input(input_file)
    msv = read_output(output_file)

    # ===== 初始图（只生成一次）=====
    plot_initial(poly1, poly2, save_dir)

    total_moves = len(moves)
    sample_size = min(1000, total_moves)

    # 随机选 1000 个
    sampled_indices = sorted(random.sample(range(total_moves), sample_size))

    # 绘制对比图
    for step, idx in enumerate(sampled_indices, start=1):
        move = moves[idx]
        sep_vec = msv[idx]
        move_idx = idx + 1  # 原始编号从1开始

        plot_compare(poly1, poly2, move, sep_vec, step, move_idx, save_dir)

print("✅ 全部完成（简洁命名版 s*_m*.png）")