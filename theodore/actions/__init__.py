from .actions import ActionFactory
from .theoinp import TheodoreInput
from .analyze_tden import AnalyzeTden, AnalyzeTdenUnr, AnalyzeTdenEs2Es
from .analyze_tden_soc import AnalyzeTdenSoc
from .analyze_sden import AnalyzeSden
from .analyze_NOs import AnalyzeNOs
from .parse_libwfa import ParseLibwfa
from .plot_vist import PlotVist
from .nics_bq import NICSBq
from .plot_OmFrag import PlotOmFrag
from .plot_Om_bars import PlotOmBars
from .plot_frag_decomp import PlotFragDecomp
from .plot_graph import PlotGraph
from .plot_graph_nx import PlotGraphNx
from .jmol_MOs import JMolMOs
from .jmol_vibs import JMolVibs
from .vmd_plots import VMDPlots
from .draw_moments import DrawMoments
from .babel import Babel
from .cc_opt import CCOpt
from .cc_check import CCCheck
from .extract_molden import ExtractMolden
from .spectrum import Spectrum
from .tden_OV import TDenOv
from .convert_table import ConvertTable
from .dgrid_prep import DGridPrep
from .fcd import FCD


from .. import theo_header


settings = {
            'description': theo_header.ret_header(),
            'logo': None,
            'error_order': ['logo', 'description', 'args', 'usage', 'space', 'comment', 'space', 'error'],
            'arg_format': {
                'name': 20,
                'comment': 50,
                'seperator': ' | ',
            },
            'subparser_args': {
                'title': 'Actions: %s',
            },
            'subparser_format': {
                'name': 20,
                'comment': 50,
                },
}

run = ActionFactory.from_commandline(description=settings, as_parser=True)
