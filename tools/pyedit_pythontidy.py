# Copyright (C) 2008-today The SG++ project
# This file is part of the SG++ project. For conditions of distribution and
# use, please see the copyright notice provided with SG++ or at 
# sgpp.sparsegrids.org


"""
This code is public domain.
The original author is Bear Huang (http://bear330.wordpress.com/).
"""
if False:
    from org.python.pydev.editor import PyEdit #@UnresolvedImport
    cmd = 'command string'
    editor = PyEdit

assert cmd is not None
assert editor is not None

if cmd == 'onCreateActions':
    # from org.eclipse.jface.action import Action
    from org.python.pydev.editor.actions import PyAction
    from org.python.pydev.core.docutils import PySelection
    from java.lang import Runnable
    from org.eclipse.swt.widgets import Display
    from org.eclipse.jface.text import IDocument
    from org.eclipse.jface.text import TextSelection

    from java.io import FileWriter
    import java.lang.Exception

    FORMAT_ACTION_DEFINITION_ID = "org.python.pydev.editor.actions.pyFormatStd"
    FORMAT_ACTION_ID = "org.python.pydev.editor.actions.navigation.pyFormatStd"

    class PythonTidyAction(PyAction):

        def __init__(self, *args, **kws):
            PyAction.__init__(self, *args, **kws)

        def run(self):
            import tempfile
            import os

            try:
                ps = PySelection(editor)
                doc = ps.getDoc()
                startLine = ps.getStartLineIndex()

                p1 = tempfile.mktemp()
                p2 = tempfile.mktemp()
                f1 = FileWriter(p1)

                formatAll = False
                if ps.getTextSelection().getLength() == 0:
                    # format all.
                    c = doc.get()
                    f1.write(c)
                    formatAll = True
                else:
                    # format selection.
                    #c = ps.getSelectedText()
                    #f1.write(ps.getSelectedText())
                    print("Format selected text is not supported yet.")
                    f1.write("")
                    # A kind of solution is to insert a special comment in
                    # front and end of selection text, pythontidy it, and
                    # extract text according that comment. 

                f1.close()
                os.system('PythonTidy.py "%s" "%s"' % (p1, p2))
                f2 = open(p2, "r")
                result = f2.read()
                f2.close()

                os.remove(p1)
                os.remove(p2)

                if startLine >= doc.getNumberOfLines():
                    startLine = doc.getNumberOfLines() - 1

                if formatAll:
                    doc.set(result)
                else:
                    #doc.replace(doc.getLineOffset(startLine), 0, result)
                    pass

                sel = TextSelection(doc, doc.getLineOffset(startLine), 0)
                self.getTextEditor().getSelectionProvider().setSelection(sel)
            except java.lang.Exception as e:
                self.beep(e)

    def bindInInterface():
        act = PythonTidyAction()

        act.setActionDefinitionId(FORMAT_ACTION_DEFINITION_ID)
        act.setId(FORMAT_ACTION_ID)
        try:
            editor.setAction(FORMAT_ACTION_ID, act)
        except:
            pass

    class RunInUi(Runnable):
        '''Helper class that implements a Runnable (just so that we
        can pass it to the Java side). It simply calls some callable.
        '''

        def __init__(self, c):
            self.callable = c
        def run(self):
            self.callable ()

    def runInUi(callable):
        '''
        @param callable: the callable that will be run in the UI
        '''
        Display.getDefault().asyncExec(RunInUi(callable))

    runInUi(bindInInterface)
