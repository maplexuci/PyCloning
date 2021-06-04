from tkinter import *
from tkinter import scrolledtext

root = Tk()

frame = Frame(root)
frame.pack(fill=BOTH, expand=True)
frame.update()

textEditor = scrolledtext.ScrolledText(frame, wrap=WORD,
                                        font=("ConsolasSeq", 12), bg="#f5feff")

textEditor.pack(expand=True, fill=BOTH, anchor=W)

seq = """AAATATCGATAGCTTTAAAGCGCCCCCTTTGGGCGGCC
CAAAAACTTCGGGCGTCTTTCAAAACGGCATATATACG
GGCGGCGGATATTAAACTTTTCGGGAGCGTCATTATCG
GAGCGAGGGCTTCTTAGGAGCTTTGGCGGGAGGATCTT
CTCTGGAGCGGGCGCTATTAACATTTCTATATTTAAAA
CTTCTCTTCTTAAAACTTTTCTCTATATATAAATTCTT
CTTTTTTCCCATATTATATCTTTCTTCTCTTTTTTCTC
TTATTATATATAAATTCTTTCTTATAACCCTTATAAAA
AAAAAACTTCTTTTAAACCTTATTC
"""

textEditor.insert('1.0', seq)

def selection(event):
    """select two lines at a time"""
    # Keep tracking where the cursor is when it is clicked
    start_position = textEditor.index(INSERT)
    print(start_position)

    # Keep tracking where the cursor is when it is dragged.


    # Coumpute which character is in the bounding box.


    # Add sel tag to the those characters.
    
    pass


# Bind <ButtonPress-1>, <B1-Motion> and <ButtonRelease-1> with selection function.
textEditor.bindtags(('Text', 'post-class-bindings', '.', 'all'))
textEditor.bind_class("post-class-bindings", "<Button-1>", selection)
# textEditor.bindtags(('Text', '.textEditor', '.', 'all'))
# textEditor.bind("<Button-1>", selection)

root.mainloop()