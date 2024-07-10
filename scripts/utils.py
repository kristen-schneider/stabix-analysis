def read_colors(colors_file):
    '''
    read comma delimited colors file and store in a dict
    :param colors_file:
    :return:
    '''
    colors = {}
    with open(colors_file, 'r') as f:
        for line in f:
            codec, color = line.strip().split(',')
            colors[codec] = color

    return colors