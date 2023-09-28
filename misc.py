def is_int( value):
    try:
        int(value)
        return True
    except ValueError:
        print(f"Error: {value} is not an integer.")
        return False

class OptionNotFoundError(Exception):
    pass

def is_option(value, options):
    try:
        if value not in options:
            raise OptionNotFoundError(f"Error: {value} is not in {options}.")
        else:
            return ValueError
    except OptionNotFoundError as e:
        print(e)
