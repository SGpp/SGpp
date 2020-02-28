# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at
# sgpp.sparsegrids.org

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

            dr, _, _ = select.select([sys.stdin], [], [], 0)
            while dr == [] and timeout > 0:
                time.sleep(1)
                timeout = timeout - 1
                sys.stdout.write("\b"*len("%02d"%(timeout+1)) + "%02d"%(timeout))
                sys.stdout.flush()
                dr, _, _ = select.select([sys.stdin], [], [], 0)

            if dr != []:
                keypressed = True
        finally:
            try:
                # switch to normal terminal
                termios.tcsetattr(fd, termios.TCSAFLUSH, old_term)
            except:
                sys.stdout.write("\nAttention: Terminal unacessable. If SGpp was built within Eclipse, everything is fine!\n")
                sys.stdout.flush()
                return keypressed

    sys.stdout.write("\n")
    sys.stdout.flush()
    return keypressed
