# Copyright (C) 2009  Dirk Pflueger (pflueged@in.tum.de)
# Version 1.0, 2009/05/16

import sys, time

def countdown(timeout):
    '''A countdown function, asking for a keypress.
    It returns whether a key was pressed during countdown or not.
    Works on both linux on windows. Timeout in seconds.'''
    sys.stdout.write("Press any key to abort... %02d" % (timeout))
    sys.stdout.flush()
    keypressed = False
    try:
        # Windows
        import msvcrt
        while not msvcrt.kbhit() or timeout == 0:
            time.sleep(1)
            timeout = timeout - 1
            sys.stdout.write("\b"*len("%02d"%(timeout+1)) + "%02d"%(timeout))
            sys.stdout.flush()
        if msvcrt.kbhit():
            keypressed = True
    except:
        # Linux
        try:
            import termios, select
            fd = sys.stdin.fileno()
            new_term = termios.tcgetattr(fd)
            old_term = termios.tcgetattr(fd)
            # new terminal setting unbuffered
            new_term[3] = (new_term[3] & ~termios.ICANON & ~termios.ECHO)
            # switch to unbuffered terminal
            termios.tcsetattr(fd, termios.TCSAFLUSH, new_term)

            dr,dw,de = select.select([sys.stdin], [], [], 0)
            while dr == [] and timeout > 0:
                time.sleep(1)
                timeout = timeout - 1
                sys.stdout.write("\b"*len("%02d"%(timeout+1)) + "%02d"%(timeout))
                sys.stdout.flush()
                dr,dw,de = select.select([sys.stdin], [], [], 0)
            if dr <> []:
                keypressed = True
        finally:
            # switch to normal terminal
            termios.tcsetattr(fd, termios.TCSAFLUSH, old_term)

    sys.stdout.write("\n")
    sys.stdout.flush()
    return keypressed
