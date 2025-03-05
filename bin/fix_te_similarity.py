#!/usr/bin/env python3
import sys

def add_batch_function(filename):
    with open(filename, 'r') as f:
        content = f.read()
    
    # 查找是否已经有batch函数
    if "def calculate_sequence_similarity_batch(" in content:
        print(f"{filename} 已经包含 calculate_sequence_similarity_batch 函数")
        return False
    
    # 找一个好的位置添加函数
    insertion_point = content.find("def main(")
    if insertion_point == -1:
        print(f"在 {filename} 中没有找到合适的插入点")
        return False
    
    # 添加批处理函数
    batch_function = """
def calculate_sequence_similarity_batch(batch_tasks):
    \"\"\"批量计算序列相似度\"\"\"
    batch_results = []
    
    for te_id, inserted_seq, original_seq, info in batch_tasks:
        # 简化的相似度计算
        min_len = min(len(original_seq), len(inserted_seq))
        
        matches = 0
        for i in range(min_len):
            if original_seq[i] == inserted_seq[i]:
                matches += 1
        
        similarity = (matches / min_len) * 100 if min_len > 0 else 0
        identity = (matches / len(original_seq)) * 100 if len(original_seq) > 0 else 0
        
        batch_results.append((te_id, {
            'similarity': similarity,
            'identity': identity,
            'original_length': len(original_seq),
            'inserted_length': len(inserted_seq),
            'original_id': info['original_id'],
            'nested_in': info.get('nested_in', 'None')
        }))
    
    return batch_results

"""
    
    # 插入函数
    new_content = content[:insertion_point] + batch_function + content[insertion_point:]
    
    # 写入修改后的文件
    with open(filename, 'w') as f:
        f.write(new_content)
    
    print(f"已向 {filename} 添加 calculate_sequence_similarity_batch 函数")
    return True

def fix_parallel_processing(filename):
    with open(filename, 'r') as f:
        content = f.read()
    
    # 查找calculate_similarity_metrics_parallel函数
    function_start = content.find("def calculate_similarity_metrics_parallel(")
    if function_start == -1:
        print(f"在 {filename} 中没有找到 calculate_similarity_metrics_parallel 函数")
        return False
    
    # 查找函数结束
    next_def = content.find("def", function_start + 1)
    if next_def == -1:
        next_def = len(content)
    
    # 提取函数内容
    function_content = content[function_start:next_def]
    
    # 创建修复后的函数
    fixed_function = """def calculate_similarity_metrics_parallel(library_seqs, inserted_seqs, num_processes=None):
    \"\"\"并行计算TE序列相似度指标\"\"\"
    if num_processes is None:
        num_processes = get_optimal_processes(task_type='cpu_bound')
    
    # 准备计算任务
    tasks = []
    
    # 预处理：收集所有需要比较的序列对
    for te_id, inserted_data in inserted_seqs.items():
        inserted_seq = inserted_data['seq']
        desc = inserted_data['desc']
        
        # 从描述中提取信息
        info = {'original_id': te_id.rsplit('_', 1)[0], 'nested_in': 'None'}
        for item in desc.split():
            if '=' in item:
                key, value = item.split('=', 1)
                info[key] = value
        
        # 获取原始TE序列
        original_id = info['original_id']
        if original_id in library_seqs:
            original_seq = library_seqs[original_id]['seq']
            
            # 如果是反向互补插入，需要调整原始序列
            if info.get('strand') == '-':
                original_seq = reverse_complement(original_seq)
            
            tasks.append((te_id, inserted_seq, original_seq, info))
    
    # 批处理任务
    total_tasks = len(tasks)
    batch_size = max(1, min(100, total_tasks // (num_processes * 2)))
    
    # 对任务进行分组
    batched_tasks = []
    for i in range(0, len(tasks), batch_size):
        batched_tasks.append(tasks[i:i+batch_size])
    
    # 使用进程池并行处理批次任务
    print(f"TE序列相似度计算：将{total_tasks}个任务分为{len(batched_tasks)}批进行处理")
    
    # 使用map而非imap_unordered来避免某些版本的Python中的序列化问题
    from parallel_processing import calculate_sequence_similarity_batch
    with mp.Pool(processes=num_processes) as pool:
        all_batch_results = pool.map(calculate_sequence_similarity_batch, batched_tasks)
        
        # 收集结果
        metrics = {}
        for batch_results in all_batch_results:
            for te_id, data in batch_results:
                metrics[te_id] = data
    
    return metrics
"""
    
    # 替换原始函数
    new_content = content[:function_start] + fixed_function + content[next_def:]
    
    # 写入修复后的文件
    with open(filename, 'w') as f:
        f.write(new_content)
    
    print(f"已修复 {filename} 中的 calculate_similarity_metrics_parallel 函数")
    return True

if __name__ == "__main__":
    # 先修复parallel_processing.py
    add_batch_function("parallel_processing.py")
    fix_parallel_processing("parallel_processing.py")
