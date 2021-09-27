"""
Copyright 2020 William Rochira at York Structural Biology Laboratory
"""

NO_CONTEXT_MESSAGE = 'NO CONTEXT'
SC_INCOMPLETE_STRING = 'INCOMPLETE SIDECHAIN'

RESOLUTION_BIN_NAMES = ('<10', '10-20', '20-30', '30-40', '40-50', '50-60', '60-70', '70-80', '80-90', '>90', 'All')

METRIC_NAMES = ('Ramachandran Score', 'Rotamer Score', 'Avg B-factor', 'Max B-factor', 'Std B-factor', 'Residue Fit', 'Mainchain Fit', 'Sidechain Fit', 'wRMSD', 'False Negative Count')
METRIC_POLARITIES = (+1, -1, -1, -1, -1, -1, -1, -1, -1, -1)
METRIC_SHORTNAMES = ('Rama', 'Rota', 'Avg B', 'Max B', 'Std B', 'Res. Fit', 'MC Fit', 'SC Fit', 'wRMSD', 'FN Count')
METRIC_DISPLAY_TYPES = ('D', 'D', 'C', 'C', 'C', 'C', 'C', 'C', 'C', 'C')

REPORT_METRIC_IDS = (0, 1, 2, 3, 6, 7)
REPORT_RESIDUE_VIEW = 'Grid'


COVARIANCE_DATA_TABLE = """ <table {}>
  <tr>
    <th>Outlier no.</th>
    <th>Residue no.</th>
    <th>wRMSD</th>
    <th>FN Count</th>
    <th>Show fix</th>
  </tr>
  {}
  </table> """

COVARIANCE_DATA_ROW="""<tr>
    <td>{}</td>
    <td>{}</td>
    <td>{}</td>
    <td>{}</td>
    <td><button id={}_{} type="button" onclick="showCovFixTable(this.id)"; class="btn btn-default">Show</button></td>
  </tr>"""


COVARIANCE_DATA_NO_OUTLIER_ROW="""<tr>
    <td colspan="5">No outliers found.</td>
  </tr>"""


COVARIANCE_FIX_TABLE = """ <table id="covariance_fix_table_{}_{}" style="display:none">
  <tr>
    <th>Current Residue</th>
    <th>New Residue</th>
  </tr>
  {}
  </table> """

COVARIANCE_FIX_ROW ="""<tr>
    <td>{}</td>
    <td>{}</td>
  </tr>"""

COVARIANCE_NO_FIX_ROW = """<tr>
    <td colspan="2">No fix found</td>
  </tr>"""
