import matplotlib.pyplot as plt

# 解析坐标
def parse_polygon(line, n):
    nums = list(map(float, line.split()))
    return [(nums[2*i], nums[2*i+1]) for i in range(n)]

# 绘制
def draw_polygon(points, color, label):
    xs = [p[0] for p in points] + [points[0][0]]
    ys = [p[1] for p in points] + [points[0][1]]
    plt.plot(xs, ys, color=color, label=label)

# 平移
def translate_polygon(points, dx, dy):
    return [(x + dx, y + dy) for x, y in points]

def main():
    print("请输入两个多边形的顶点数（例如：17 15）：")
    n1, n2 = map(int, input().split())

    print("请输入第一个多边形坐标：")
    poly1 = parse_polygon(input(), n1)

    print("请输入第二个多边形坐标：")
    poly2 = parse_polygon(input(), n2)

    # ⭐ 当前状态（关键！）
    current_poly2 = poly2.copy()

    while True:
        print("\n请输入移动向量 dx dy（输入 q 退出）：")
        s = input()
        if s.lower() == 'q':
            break

        dx, dy = map(float, s.split())

        # ⭐ 在“当前状态”基础上移动
        current_poly2 = translate_polygon(current_poly2, dx, dy)

        plt.clf()
        draw_polygon(poly1, 'blue', 'Polygon 1 (固定)')
        draw_polygon(current_poly2, 'red', 'Polygon 2 (累计移动)')

        plt.legend()
        plt.title(f"累计移动: dx={dx}, dy={dy}")
        plt.axis('equal')
        plt.grid(True)
        plt.pause(0.1)

    plt.show()

if __name__ == "__main__":
    main()