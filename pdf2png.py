from pdf2image import convert_from_bytes, convert_from_path
import os
def pdf2png(pdf_file):
    pdf_filename=pdf_file.split("/")[-1]

    if ".pdf" == pdf_file[-4]:
        png_dir="%s_png"%(pdf_file[:-4])
    else:
        raise Exception

    os.makedirs(png_dir)
    images = convert_from_path(pdf_file)

    for i, im in enumerate(images):
        im.save("%s/%s_%s.png" % (png_dir,pdf_filename[:-4],i), "PNG")


if __name__ == '__main__':
    pdf2png("")