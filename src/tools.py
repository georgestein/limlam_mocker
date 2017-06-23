import datetime

class empty_table():
    """ 
    brief Class describing a table.    
    """
    def __init__(self):
        pass

    def copy(self):
        """ 
        @brief Creates a copy of the table.          
        """
        return copy.copy(self)

def write_time(string_in):
    fmt       = '%H:%M:%S on %m/%d/%Y'
    timestamp = datetime.datetime.now().strftime(fmt)
    bar = 72*'-'
    print '\n\n'+bar
    print string_in
    print 'Time:      '+timestamp
    print bar+'\n'

    return
