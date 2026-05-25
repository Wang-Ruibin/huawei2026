def compare_files(file1, file2):
    with open(file1, 'rb') as f1, open(file2, 'rb') as f2:
        while True:
            b1 = f1.read(4096)
            b2 = f2.read(4096)
            if b1 != b2:
                return False
            if not b1:  # 读完了
                return True

file1 = "practice_7.out"
file2 = "practice_7 copy.out"

if compare_files(file1, file2):
    print("两个文件内容完全相同！")
else:
    print("两个文件内容不同！")