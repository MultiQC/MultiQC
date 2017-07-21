from multiqc.modules.base_module import BaseMultiqcModule
from .relatedness import Relatedness2Mixin


class MultiqcModule(BaseMultiqcModule, Relatedness2Mixin):
    def __init__(self):
        super(MultiqcModule, self).__init__(
            name='VCFTools',
            anchor='vcftools',
            href='https://vcftools.github.io',
            info='is a program for working with and reporting on VCF files'
        )

        self.parse_relatedness2()
