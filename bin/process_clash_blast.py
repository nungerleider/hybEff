import sys

inpath = sys.argv[1]
outpath = inpath + '.best_mir_only.tsv'

with open(inpath) as infile:
    d = {}
    for line in infile:
        if 'miR' not in line and 'let' not in line:
            continue
        spl = line.split('\t')
        line_id = spl[0]
        score = float(spl[2])
        length = int(spl[3])
        start = int(spl[6])
        end = int(spl[7])
        read_start = int(spl[8])

        if line_id not in d:
            d[line_id] = line

        elif abs(int(d[line_id].split('\t')[6]) - start) > 16:
            if f'b{line_id}' not in d: 
                d[f'b{line_id}'] = line

            elif float(d[f'b{line_id}'].split('\t')[2]) <= score and int(d[f'b{line_id}'].split('\t')[3]) <= length:
                d[f'b{line_id}'] = line
            
            elif abs(int(d[f'b{line_id}'].split('\t')[6]) - start) > 16:
                d[f'c{line_id}'] = line
            
            elif float(d[f'b{line_id}'].split('\t')[2]) == score and int(d[f'b{line_id}'].split('\t')[3]) == length and int(d[f'b{line_id}'].split('\t')[8]) > read_start:
                d[f'b{line_id}'] = line

        elif float(d[line_id].split('\t')[2]) == score and int(d[line_id].split('\t')[3]) == length:
            if int(d[line_id].split('\t')[8]) > read_start:
                d[line_id] = line

        elif float(d[line_id].split('\t')[2]) <= score and int(d[line_id].split('\t')[3]) <= length:          
            d[line_id] = line

with open(outpath, 'w') as outfile:
    for val in d.values():
        spl = val.split('\t')
        counts = spl[0].split('_')[1]
        outfile.write(f'{counts}\t{spl[1]}\n')