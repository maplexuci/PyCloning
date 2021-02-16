from Bio.Seq import Seq
from tkinter import *
from tkinter import scrolledtext
from PIL import ImageTk, Image


def main():
    """The main app function"""
    root = Tk()
    root_window = Root(root)
    return None


class Root:
    """Root window class including the widgets and the functions for 
    responding the widgets.
    """
    def __init__(self, root):
        # Main root window configration
        self.root = root
        self.root.title("PyCloning")
        self.root.geometry("700x400")
        self.root.resizable(width=False, height=False)
        self.root.configure(bg="#ebf5fc")
        self.root.columnconfigure(0, weight=2)
        self.root.columnconfigure((1, 2), weight=1)

        # Create menu bar
        menubar = Menu(self.root)
        self.root.config(menu=menubar)

        # Create 'File' menu
        file_menu = Menu(menubar, tearoff=False)
        menubar.add_cascade(label="File", menu=file_menu)
        file_menu.add_command(label="New DNA File...", command=self.new_DNA)
        file_menu.add_command(label="New Protein File...")
        file_menu.add_separator()
        file_menu.add_command(label="Open Files...")
        file_menu.add_command(label="Open Recent File")
        file_menu.add_command(label="Close", command=quit)
        file_menu.add_separator()
        file_menu.add_command(label="Import")
        file_menu.add_separator()
        file_menu.add_command(label="Exit", command=quit)

        # Create 'Cloning' menu
        cloning_menu = Menu(menubar, tearoff=False)
        menubar.add_cascade(label="Cloning", menu=cloning_menu)
        cloning_menu.add_command(label="PCR...")
        cloning_menu.add_command(label="Overlapping PCR")
        cloning_menu.add_command(label="Mutagenesis...")
        cloning_menu.add_separator()
        cloning_menu.add_command(label="Restriction Enzyme Cloning")
        cloning_menu.add_command(label="Linear Ligation")
        cloning_menu.add_separator()
        cloning_menu.add_command(label="TA or GC Cloning")
        cloning_menu.add_command(label="TOPO® Cloning")
        cloning_menu.add_command(label="Gateway® Cloning")
        cloning_menu.add_command(label="Gibson Assembly®")
        cloning_menu.add_command(label="In-Fusion® Cloning")
        cloning_menu.add_command(label="NEBuilder® HiFi DNA Assembly")

        # Create 'Reference' menu
        ref_menu = Menu(menubar, tearoff=False)
        menubar.add_cascade(label="Reference", menu=ref_menu)
        ref_menu.add_command(label="Restriction Enzymes")
        ref_menu.add_command(label="Enzyme Database")
        ref_menu.add_separator()
        ref_menu.add_command(label="Common Features")
        ref_menu.add_command(label="Letter Codes")
        ref_menu.add_separator()
        ref_menu.add_command(label="Genetic Code Tables")
        ref_menu.add_command(label="Codon Usage Tables")

        VERSION = "V0.1"
        ver_title = Label(self.root, text="PyCloning " + VERSION, font=("Courier", 24, "bold"), bg="#ebf5fc")
        ver_title.grid(row=0, column=0, columnspan=3, padx=220, pady=10, sticky=E+W)

        canvas = Canvas(self.root, bg="#ebf5fc", highlightthickness=0, width=260, height=260)
        canvas.grid(row=1, column=0, rowspan=5, padx=3, pady=9)
        logo = ImageTk.PhotoImage(Image.open("imgs/logo.png"))
        canvas.create_image(130, 130, image=logo)

        status = Label(self.root, bd=1, relief=SUNKEN, height=3)
        status.grid(row=6, column=0, columnspan=3, sticky=W+E)

        pixelVirtual = PhotoImage(width=1, height=1)
        btn_ndf = Button(self.root, image=pixelVirtual, width=35, height=35)
        btn_npf = Button(self.root, image=pixelVirtual, width=35, height=35)
        btn_open = Button(self.root, image=pixelVirtual, width=35, height=35)
        btn_orf = Button(self.root, image=pixelVirtual, width=35, height=35)
        btn_import = Button(self.root, image=pixelVirtual, width=35, height=35)

        lbl_ndf = Label(self.root, text="New DNA File...", bg="#ebf5fc")
        lbl_npf = Label(self.root, text="New Protein File...", bg="#ebf5fc")
        lbl_open = Label(self.root, text="Open", bg="#ebf5fc")
        lbl_orf = Label(self.root, text="Open Recent File", bg="#ebf5fc")
        lbl_import = Label(self.root, text="Import", bg="#ebf5fc")

        btn_ndf.grid(row=1, column=1, pady=5)
        btn_npf.grid(row=2, column=1, pady=5)
        btn_open.grid(row=3, column=1, pady=5)
        btn_orf.grid(row=4, column=1, pady=5)
        btn_import.grid(row=5, column=1, pady=5)

        lbl_ndf.grid(row=1, column=2, pady=10, sticky=W)
        lbl_npf.grid(row=2, column=2, pady=10, sticky=W)
        lbl_open.grid(row=3, column=2, pady=10, sticky=W)
        lbl_orf.grid(row=4, column=2, pady=10, sticky=W)
        lbl_import.grid(row=5, column=2, pady=10, sticky=W)

        self.root.mainloop()

    def hide(self):
        """Hide the root window."""
        self.root.withdraw()

    def show(self):
        """Show the root window from the hide status"""
        self.root.update()
        self.root.deiconify()

    def onClosing(self, window):
        """Respond to the toplevle window closing event:
        Close the current window and show the root window.
        """
        window.destroy()
        self.show()

    def new_DNA(self):
        """Create new toplevel window for new DNA file

        Returns:
            [type]: [description]
        """
        self.hide()
        self.dna_window = Toplevel()
        self.dna_window.title("New DNA File")
        self.dna_window.geometry("500x400")
        # 'grab_set()' and `focus()` ensure the toplevel is active (focused).
        self.dna_window.grab_set()
        self.dna_window.focus()

        lbl = Label(self.dna_window, text="Input your sequence:")
        lbl.pack(padx=10, pady=(10, 0), anchor=W)

        # Create a scrolledtext widget.
        self.new_dna_seq = scrolledtext.ScrolledText(self.dna_window, wrap=WORD, width=80, height=15, font=("Courier New", 11))
        # 'expand=True' and 'fill=BOTH' ensure that
        # the Text widget change size along with window resizing.
        self.new_dna_seq.pack(padx=10, expand=True, fill=BOTH, anchor=W)
        self.new_dna_seq.focus_set()
        self.new_dna_seq.bind('<<Modified>>', self.changed)

        # Add another Label widget to display the sequence length.
        self.seq_len_lbl = Label(self.dna_window)
        self.seq_len_lbl.pack(padx=10, anchor=W)

        # Add Label and Entry widget in a Frame for entering filename.
        self.filename_Frame = Frame(self.dna_window)
        Label(self.filename_Frame, text="File Name: ").pack(side=LEFT)
        self.name_var = StringVar()
        self.seq_filename = Entry(self.filename_Frame, textvariable=self.name_var)
        self.seq_filename.pack(side=LEFT)
        self.filename_Frame.pack(padx=10, anchor=W)

        # Respond to the 'Cancel' button.
        btn_cancel = Button(self.dna_window, text="Cancel", width=10,
                            command=lambda: self.onClosing(self.dna_window))
        btn_cancel.pack(padx=10, pady=10, side=RIGHT)

        # Add 'OK' button to read sequence
        btn_ok = Button(self.dna_window, text="OK", width=10,
                            command=self.readSeq)
        btn_ok.pack(padx=10, pady=10, side=RIGHT)

        #  Respond to 'close window' event
        self.dna_window.protocol("WM_DELETE_WINDOW",
                            lambda: self.onClosing(self.dna_window))

        return None

    def changed(self, seq=None):
        """Function to keep track the changes in Text and reflect the changes in the Label"""
        text = ''
        self.flag = self.new_dna_seq.edit_modified()
        if self.flag == 1:     # prevent from getting called twice
            text = str(len(self.new_dna_seq.get(1.0, END))) + " bp"
            self.seq_len_lbl.config(text=text)
            self.flag = self.new_dna_seq.edit_modified(False)  # Reset the flag to 0.

    def get_EntryContent(self, text_var):
        """Get the content in an Entry as string.

        Args:
            name_var ([type]): The variable type defined for the Entry.

        Returns:
            Str: The content in an Entry widget.
        """
        self.text_var = text_var
        return self.text_var.get()

    def readSeq(self):
        self.dnaSeq = Seq(self.new_dna_seq.get(1.0, END))  # '1.0' means the first row (number), the first column (index).
        self.dna_window.destroy()
        self.seq_window = Toplevel()
        self.seq_window.title(self.get_EntryContent(self.name_var))
        self.seq_window.state('zoomed')  # Make the window maximized.
        self.seq_window.protocol("WM_DELETE_WINDOW",
                            lambda: self.onClosing(self.seq_window))

        # Create a scrolledtext widget.
        self.dna_seq = scrolledtext.ScrolledText(self.seq_window, wrap=WORD, font=("Courier New", 11))
        # 'expand=True' and 'fill=BOTH' ensure that
        # the Text widget change size along with window resizing.
        self.dna_seq.pack(padx=10, expand=True, fill=BOTH, anchor=W)
        self.dna_seq.insert(1.0, self.dnaSeq)
        # self.dna_seq.bind('<<Modified>>', self.changed)

        


main()
