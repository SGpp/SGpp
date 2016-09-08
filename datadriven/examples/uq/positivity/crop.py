import subprocess
import os
import shutil

# files that need cropping are added to filesToCrop list
for pdffile in os.listdir("figures"):
    pdffile = os.path.join("figures", pdffile)
    if pdffile.endswith(".pdf"):
        subprocess.call(["pdfcrop", pdffile, pdffile])
