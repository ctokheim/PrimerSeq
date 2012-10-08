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

