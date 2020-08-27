FROM rocker/r-ver:3.6.3

RUN mkdir /home/MultiGWAS

COPY INSTALL_Rlibs.R /home/MultiGWAS

CMD R -e "source ('/home/MultiGWAS/INSTALL_Rlibs.R')"
