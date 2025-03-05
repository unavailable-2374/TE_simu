#!/usr/bin/env python3
import sys

def fix_file(filename):
    with open(filename, 'r') as f:
        content = f.read()
    
    # 查找calculate_gc_content函数
    function_start = content.find("def calculate_gc_content(sequence")
    if function_start == -1:
        print(f"在 {filename} 中没有找到 calculate_gc_content 函数")
        return False
    
    # 查找函数结束
    next_def = content.find("def", function_start + 1)
    if next_def == -1:
        next_def = len(content)
    
    # 提取函数内容
    function_content = content[function_start:next_def]
    
    # 检查函数是否已经修复
    if "isinstance(sequence, (bytes, bytearray))" in function_content:
        print(f"{filename} 中的 calculate_gc_content 函数已经修复")
        return False
    
    # 创建修复后的函数
    fixed_function = """def calculate_gc_content(sequence):
    \"\"\"计算序列的GC含量，支持字符串和字节对象\"\"\"
    if not sequence:
        return 0
    
    # 检查序列类型并相应处理
    if isinstance(sequence, (bytes, bytearray)):
        g_count = sequence.count(ord('G')) + sequence.count(ord('g'))
        c_count = sequence.count(ord('C')) + sequence.count(ord('c'))
        return (g_count + c_count) / len(sequence) * 100
    else:
        # 原始字符串处理方式
        g_count = sequence.upper().count('G')
        c_count = sequence.upper().count('C')
        return (g_count + c_count) / len(sequence) * 100
"""
    
    # 替换原始函数
    new_content = content[:function_start] + fixed_function + content[next_def:]
    
    # 写入修复后的文件
    with open(filename, 'w') as f:
        f.write(new_content)
    
    print(f"已修复 {filename} 中的 calculate_gc_content 函数")
    return True

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("用法: python fix_calculate_gc_content.py <filename>")
        sys.exit(1)
    
    fix_file(sys.argv[1])
