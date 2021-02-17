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
    """Root window class including the widgets and the functions
    for responding the widgets.
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

        # Create menubar first
        self.menubar = Menu(self.root)
        self.root.config(menu=self.menubar)

        # Create 'File' menu inside the menubar
        # 'tearoff=False' ensures the menu is not a floatting menu
        self.file_menu = Menu(self.menubar, tearoff=False)
        # Adding cascade to 'File' menu.
        self.menubar.add_cascade(label="File", menu=self.file_menu)
        self.file_menu.add_command(label="New DNA File...", command=self.new_DNA)
        self.file_menu.add_command(label="New Protein File...")
        self.file_menu.add_separator()
        self.file_menu.add_command(label="Open Files...")
        self.file_menu.add_command(label="Open Recent File")
        self.file_menu.add_command(label="Close", command=quit)
        self.file_menu.add_separator()
        self.file_menu.add_command(label="Import")
        self.file_menu.add_separator()
        self.file_menu.add_command(label="Exit", command=quit)

        # Create 'Cloning' menu
        self.cloning_menu = Menu(self.menubar, tearoff=False)
        self.menubar.add_cascade(label="Cloning", menu=self.cloning_menu)
        self.cloning_menu.add_command(label="PCR...")
        self.cloning_menu.add_command(label="Overlapping PCR")
        self.cloning_menu.add_command(label="Mutagenesis...")
        self.cloning_menu.add_separator()
        self.cloning_menu.add_command(label="Restriction Enzyme Cloning")
        self.cloning_menu.add_command(label="Linear Ligation")
        self.cloning_menu.add_separator()
        self.cloning_menu.add_command(label="TA or GC Cloning")
        self.cloning_menu.add_command(label="TOPO® Cloning")
        self.cloning_menu.add_command(label="Gateway® Cloning")
        self.cloning_menu.add_command(label="Gibson Assembly®")
        self.cloning_menu.add_command(label="In-Fusion® Cloning")
        self.cloning_menu.add_command(label="NEBuilder® HiFi DNA Assembly")

        # Create 'Reference' menu
        self.ref_menu = Menu(self.menubar, tearoff=False)
        self.menubar.add_cascade(label="Reference", menu=self.ref_menu)
        self.ref_menu.add_command(label="Restriction Enzymes")
        self.ref_menu.add_command(label="Enzyme Database")
        self.ref_menu.add_separator()
        self.ref_menu.add_command(label="Common Features")
        self.ref_menu.add_command(label="Letter Codes")
        self.ref_menu.add_separator()
        self.ref_menu.add_command(label="Genetic Code Tables")
        self.ref_menu.add_command(label="Codon Usage Tables")

        self.VERSION = "V0.1"
        self.ver_title = Label(self.root, text="PyCloning " + self.VERSION, font=("Courier", 24, "bold"), bg="#ebf5fc")
        self.ver_title.grid(row=0, column=0, columnspan=3, padx=220, pady=10, sticky=E+W)

        self.canvas = Canvas(self.root, bg="#ebf5fc", highlightthickness=0, width=260, height=260)
        self.canvas.grid(row=1, column=0, rowspan=5, padx=3, pady=9)
        self.logo = ImageTk.PhotoImage(Image.open("imgs/logo.png"))
        self.canvas.create_image(130, 130, image=self.logo)

        self.status = Label(self.root, bd=1, relief=SUNKEN, height=3)
        self.status.grid(row=6, column=0, columnspan=3, sticky=W+E)

        self.pixelVirtual = PhotoImage(width=1, height=1)
        self.btn_ndf = Button(self.root, image=self.pixelVirtual, width=35, height=35)
        self.btn_npf = Button(self.root, image=self.pixelVirtual, width=35, height=35)
        self.btn_open = Button(self.root, image=self.pixelVirtual, width=35, height=35)
        self.btn_orf = Button(self.root, image=self.pixelVirtual, width=35, height=35)
        self.btn_import = Button(self.root, image=self.pixelVirtual, width=35, height=35)

        self.lbl_ndf = Label(self.root, text="New DNA File...", bg="#ebf5fc")
        self.lbl_npf = Label(self.root, text="New Protein File...", bg="#ebf5fc")
        self.lbl_open = Label(self.root, text="Open", bg="#ebf5fc")
        self.lbl_orf = Label(self.root, text="Open Recent File", bg="#ebf5fc")
        self.lbl_import = Label(self.root, text="Import", bg="#ebf5fc")

        self.btn_ndf.grid(row=1, column=1, pady=5)
        self.btn_npf.grid(row=2, column=1, pady=5)
        self.btn_open.grid(row=3, column=1, pady=5)
        self.btn_orf.grid(row=4, column=1, pady=5)
        self.btn_import.grid(row=5, column=1, pady=5)

        self.lbl_ndf.grid(row=1, column=2, pady=10, sticky=W)
        self.lbl_npf.grid(row=2, column=2, pady=10, sticky=W)
        self.lbl_open.grid(row=3, column=2, pady=10, sticky=W)
        self.lbl_orf.grid(row=4, column=2, pady=10, sticky=W)
        self.lbl_import.grid(row=5, column=2, pady=10, sticky=W)

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
        self.new_dna_seq = scrolledtext.ScrolledText(
                                self.dna_window, wrap=WORD,
                                width=80, height=15, font=("Courier New", 11)
                                )
        # 'expand=True' and 'fill=BOTH' ensure that
        # the Text widget change size along with window resizing.
        self.new_dna_seq.pack(padx=10, expand=True, fill=BOTH, anchor=W)
        self.new_dna_seq.focus_set()
        # Bind the scrooedtext widget with 'self.changed' function
        # to keep track if the content in the widget is 'Modified'.
        self.new_dna_seq.bind('<<Modified>>', self.changed)

        # Add another Label widget to display the sequence length.
        self.seq_len_lbl = Label(self.dna_window)
        self.seq_len_lbl.pack(padx=10, anchor=W)

        # Add Label and Entry widget in a Frame for entering filename.
        # First, create Frame widget
        filename_Frame = Frame(self.dna_window)
        # Create text Label
        Label(filename_Frame, text="File Name: ").pack(side=LEFT)
        # Define and use variable type for Entry widget.
        self.name_var = StringVar()
        seq_filename = Entry(filename_Frame, textvariable=self.name_var)
        seq_filename.pack(side=LEFT)
        # Pack the Frame in the end after all widgets in Frame are packed.
        filename_Frame.pack(padx=10, anchor=W)

        # Respond to the 'Cancel' button.
        btn_cancel = Button(self.dna_window, text="Cancel", width=10,
                            command=lambda: self.onClosing(self.dna_window))
        btn_cancel.pack(padx=10, pady=10, side=RIGHT)
        # Add 'OK' button to read sequence
        btn_ok = Button(self.dna_window, text="OK", width=10, state=DISABLED,
                        command=self.readSeq)
        btn_ok.pack(padx=10, pady=10, side=RIGHT)
        #  Respond to 'close window' event
        self.dna_window.protocol("WM_DELETE_WINDOW",
                                 lambda: self.onClosing(self.dna_window))

        return None

    def changed(self, event):
        """Function to keep track the changes in Text and reflect the changes in the Label"""
        # Only display the length of sequence when it is not empty.
        text = ''
        text_len = len(self.new_dna_seq.get(1.0, END).rstrip())
        if text_len != 0:
            text = f"{str(text_len)} bp"
        else:
            text = ''

        # Display the sequence length on the Label.
        self.seq_len_lbl.config(text=text)
        
        # Reset the modified flag to False, for detecting new modification.
        self.new_dna_seq.edit_modified(False)

    def get_EntryContent(self, text_var):
        """Get the content in an Entry as string.

        Args:
            name_var (var): The variable type defined for the Entry.

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
        self.dna_seq = scrolledtext.ScrolledText(self.seq_window, wrap=WORD,
                                                 font=("Courier New", 11))
        # 'expand=True' and 'fill=BOTH' ensure that
        # the Text widget change size along with window resizing.
        self.dna_seq.pack(padx=10, expand=True, fill=BOTH, anchor=W)
        self.dna_seq.insert(1.0, self.dnaSeq)


if __name__ == '__main__':
    main()
