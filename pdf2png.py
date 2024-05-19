from pdf2image import convert_from_bytes, convert_from_path
import os
import re
import logging

logging.basicConfig(
    level=logging.INFO,
    format='[%(asctime)s] [%(levelname)s] [%(name)s]  %(funcName)s  %(lineno)d- %(message)s')

logger = logging.getLogger(__name__)
def pdf2png(pdf_file):
    pdf_filename=pdf_file.split("/")[-1]

    if pdf_file.endswith(".pdf"):
        png_dir="%s_png"%(pdf_file[:-4])
    else:
        logger.error("pdf file name must end with .pdf")
        raise Exception
    if os.path.exists(png_dir):
        if os.path.exists("%s/OK"%(png_dir)):
            logger.info("%s png already exists"%(png_dir))
            return
    os.makedirs(png_dir, exist_ok=True)
    images = convert_from_path(pdf_file)

    for i, im in enumerate(images):
        im.save("%s/%s_%s.png" % (png_dir,pdf_filename[:-4],i), "PNG")
    with open("%s/OK"%(png_dir),"w") as fp:
        fp.write("OK")


def interfacial_reaction_pdf2png(outputs_dir):
    for file in os.listdir(outputs_dir):
        if file=="reaction_tally.pdf":
            logger.info("reaction_tally.pdf")
            pdf2png("%s/%s"%(outputs_dir,file))
        if file=="sink_report.pdf":
            logger.info("sink_report.pdf")
            pdf2png("%s/%s"%(outputs_dir,file))
        # if file=="species_report.pdf":#太大了，内存溢出
        #     pdf2png("%s/%s"%(outputs_dir,file))

        pathways_match=re.match(r"observed_molecule_\d+_pathways.pdf",file)
        if pathways_match:
            logger.info(pathways_match.group())
            pdf2png("%s/%s"%(outputs_dir,file)) #observed_molecule_8397_consumption_report.pdf
        consumption_match = re.match(r"observed_molecule_\d+_consumption_report.pdf", file)
        if consumption_match:
            logger.info(consumption_match.group())
            pdf2png("%s/%s" % (outputs_dir, file))



if __name__ == "__main__":
    import sys

    logger.info(sys.argv[1])
    interfacial_reaction_pdf2png(sys.argv[1])