def get_start_pos(coordinate):
    """
    get start from 'chr:start-stop'
    """
    return int(coordinate.split(":")[1].split("-")[0])


def get_end_pos(coordinate):
    """
    get stop from 'chr:start-stop'
    """
    return int(coordinate.split(":")[1].split("-")[1])


def get_pos(coordinate):
    """
    get (start, stop) from 'chr:start-stop'
    """
    return tuple(map(int, coordinate.split(":")[1].split("-")))


def get_chr(coordinate):
    """
    get chr from 'chr:start-stop'
    """
    return coordinate.split(":")[0]


def merge_list_of_dicts(list_of_dicts):
    '''
    This function mereges multiple dicts contained read counts from SAM/BAM file
    into one dictionary.
    '''
    merged_dict = {}
    for tmp_dict in list_of_dicts:
        all_keys = set(merged_dict) | set(tmp_dict)
        for key in all_keys:
            merged_dict[key] = merged_dict.get(key, 0) + tmp_dict.get(key, 0)
    return merged_dict

