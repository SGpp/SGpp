# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org

# INPUT_FILTER script for Doxygen, using default Doxyfile

# import what needed
import sys,re,os

# read in file
f = open(sys.argv[1],'r')
txt = f.read()
f.close()

# replace $HEAD$ (and, if enabled, $CURR$
if re.search('\$HEAD\$', txt):
    # read in version info
    f = os.popen('svn stat -v '+sys.argv[1])
    data = f.read()
    f.close()
    g = re.search('(\S*?)\s*?(\S*?)\s*?(\S*?)\s*?(\S*)$', data).groups();
    r_HEAD = 'Rev. %s (Last mod. in rev. %s by %s)'%(g[0],g[1],g[2]);
    r_CURR = 'Rev. '+g[0]
    txt = re.sub('version\s*\$HEAD\$', 'version '+r_HEAD, txt)
    txt = re.sub('version\s*\$CURR\$', 'version '+r_CURR, txt)

#replace $SVN_LOG$
if re.search('\$SVN_LOG\$', txt):
    # read in version info
    f = os.popen('svn stat -v '+sys.argv[1])
    data = f.read()
    f.close()
    g = re.search('(\S*?)\s*?(\S*?)\s*?(\S*?)\s*?(\S*)$', data).groups();
    try:
        rev = int(g[0])
        rev_prev = rev-4
        svn_log = 'svn log -v -r %d:%d' % (rev,rev_prev)
    except:
        svn_log = 'svn log -v -r HEAD:1 --limit 5'
    
    # read in log entries
    f = os.popen(svn_log)
    data = f.readlines()
    f.close()

    s = ""
    counter = 0
    while (counter < len(data)-1):
        if re.match('------------------------------------------------------------------------', data[counter]):
            counter += 1
            rev_info_line = data[counter].strip()
            rev_info_with_files = [rev_info_line]
            counter += 1
            while data[counter].strip() != "":
                rev_info_with_files.append(data[counter].strip())
                counter += 1
            rev_info_with_files_string = "<br>".join(rev_info_with_files)
            s += '''<tr><td colspan="2" onmouseover="this.innerHTML = '<tt>%s</tt>'" onmouseout="this.innerHTML = '<tt>%s</tt>'"><tt>%s</tt></td></tr><tr><td></td><td>''' % (
                rev_info_with_files_string, rev_info_line, rev_info_line)
            counter += 1
            while (not re.match('------------------------------------------------------------------------', data[counter])):
                if (not data[counter].isspace()):
                    s += data[counter].strip()+"<br>"
                counter += 1
            s += "</td></tr>"
## Old version without fancy onmouseover, which shows the changed files
##     s = ""
##     counter = 0
##     while (counter < len(data)-1):
##         if re.match('-', data[counter]):
##             counter += 1
##             s += '<tr><td colspan="2"><tt>%s</tt></td></tr><tr><td></td><td>' % (data[counter].strip())
##             counter += 1
##             while (not re.match('-', data[counter])):
##                 if (not data[counter].isspace()):
##                     s += data[counter].strip()+"<br>"
##                 counter += 1
##             s += "</td></tr>"
    s = '<table border="0" width="800"><tr><td width="100"></td><td width="700"></td></tr>'+s+'</table>'
            
    txt = re.sub('\$SVN_LOG\$', s, txt)

print(txt)