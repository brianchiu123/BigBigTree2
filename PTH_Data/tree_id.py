import re

def remove_species_from_file(input_filename, output_filename):
    with open(input_filename, 'r') as infile:
        tree_str = infile.read()
        # 正則表達式匹配所有的下劃線後面的英文
        cleaned_str = re.sub(r'_(\w+):', r':', tree_str)
    
    with open(output_filename, 'w') as outfile:
        outfile.write(cleaned_str)

# 使用範例
input_file = "final.ph"
output_file = "output.ph"
remove_species_from_file(input_file, output_file)
