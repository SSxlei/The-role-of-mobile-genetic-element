import os

folder_path = "."  # 替换成实际的文件夹路径

for filename in os.listdir(folder_path):
    if filename.startswith('.'):
        continue  # 跳过隐藏文件和上级目录

    new_name = filename.replace('fa', '').replace('region00', '')
    old_path = os.path.join(folder_path, filename)
    new_path = os.path.join(folder_path, new_name)

    os.rename(old_path, new_path)

print("File names have been updated successfully.")
