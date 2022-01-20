FROM python:3.6.14-stretch

ENV CNV2VCF_GIT_URL https://github.com/emedgene/cnv2vcf.git
ENV CNV2VCF_COMMIT  55ef4c8
ENV CNV2VCF_DIR     /opt/cnv2vcf

RUN cd /opt \
 && git clone $CNV2VCF_GIT_URL \
 && cd cnv2vcf/ \
 && git checkout $CNV2VCF_COMMIT \
 && pip install -r requirements.txt

RUN cd $CNV2VCF_DIR \
 && sed -i '1s/^/#!\/usr\/bin\/env python3\n\n/' bg_cnv2vcf.py \
 && chmod 777 bg_cnv2vcf.py

RUN cd $CNV2VCF_DIR \
 && sed -i '1s/^/#!\/usr\/bin\/env python3\n\n/' CONTRA_cnv2vcf.py \
 && chmod 777 CONTRA_cnv2vcf.py

RUN cd $CNV2VCF_DIR/exomeDepth \
 && sed -i '1s/^/#!\/usr\/bin\/env python3\n\n/' exomedepth_to_vcf.py \
 && chmod 777 exomedepth_to_vcf.py

ENV PATH "$CNV2VCF_DIR/exomeDepth:$PATH"
ENV PATH "$CNV2VCF_DIR:$PATH"
