import os
import re
from glob import glob
from pathlib import Path
import argparse
from decimal import Decimal, getcontext, InvalidOperation, ROUND_HALF_UP
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.lines import Line2D
from matplotlib.ticker import AutoMinorLocator
from matplotlib.transforms import blended_transform_factory

# ==========================================
# CONFIGURATION & SETTINGS (Static parameters)
# ==========================================

settings = {
	'start_dir': '.', 
	'pawley_ident': '_pawley_01_', 
	'pawley_ext': ('.txt', '.xy'),  # data file extensions to scan (any iterable of suffixes)
	'2Th_Ip_ident': '2Th_Ip', 
	'X_Yobs_ident': 'X_Yobs', 
	'Out_X_Ycalc_ident': 'Out_X_Ycalc', 
	'X_Difference_ident': 'X_Difference', 
	'figsize': (6, 4), 
	'X_Yobs_color': 'k', 
	'Out_X_Ycalc_color': 'r', 
	'2Th_Ip_colors': ['b', 'orange', 'pink', 'gold', 'red'], 
	'X_Difference_color': 'g', 
	'X_Yobs_label': 'Observed', 
	'Out_X_Ycalc_label': 'Calculated', 
	'2Th_Ip_label': 'Reflections', 
	'X_Difference_label': 'Difference', 
	'X_Yobs_markersize': 5, 
	'Out_X_Ycalc_markersize': 0, 
	'2Th_Ip_markersize': 10, 
	'X_Difference_markersize': 0, 
	'y_off_dashes': 1.5, 
	'x_lims': 'auto', 
	'legend_fontsize': 8, 
	'dpi': 300, 
	'transparent': True, 
	'extension': 'svg', 
	'title': None, 
	'title_font_weight': 'bold', 
	'title_font_size': 12, 
	'show_info': True, 
	'use_sg_from_outfile': True, 
	'use_sg_format': 'HERMANN-MAUGUIN+SUBSTANCE', 
	
	# VERTICAL LAYOUT CONTROLS (all values are axes-coordinate fractions, 0–1)
	'bragg_spacing': 0.01,          # Gap between multiple Bragg tick rows
	'tick_top_clearance': 0.02,     # Gap between data floor and the topmost Bragg tick row
	'tick_bottom_clearance': 0.02,  # Gap between the bottommost Bragg tick row and the top edge of the difference band
	'difference_band_height': 0.12, # Height of the band reserved for the difference curve (not the gap above it)
	'y_tol_top': 0.05,              # Headroom above the data ceiling (must be ≥ 0)
	'y_tol_bottom': 0.03,           # Headroom below the difference band (must be ≥ 0;
	                                # negative values push the band off the axes)
	'min_data_fraction': 0.25,      # Minimum axes fraction guaranteed to the data region
	'diff_shift_factor': 0.8,       # Max fraction of the difference band the curve may occupy.
	                                # If the natural amplitude would exceed this, the band auto-expands
	                                # so the curve fits with (1 - diff_shift_factor) of visual padding.
	'box_y_top': 0.67,              # Fallback ceiling for the first box when no legend is present
	'box_x': 0.98,                  # Fallback right-edge x-position (axes coords) when no legend
	'box_fontsize': 7,              # Default / largest font size (pt) for unit-cell info box text
	'box_fontsize_min': 5,          # Smallest font (pt) the layout will shrink to before accepting overlap
	'box_pad': 0.3,                 # Padding (in fontsize units) inside the rounded info box
	'box_gap': 0.01,                # Vertical gap (axes coords) between vertically stacked boxes
	'box_col_gap': 0.012,           # Horizontal gap (axes coords) between side-by-side boxes
	'box_legend_gap': 0.015,        # Whitespace tolerance between the legend's bottom and the first box
	'box_data_clearance': 0.02,     # Minimum clearance (axes coords) to keep above the data envelope
	'pastel_weight': 0.78,          # Blend factor toward white for info-box background tinting
	'multiply_label_y': 0.98,       # Axes-coord y for the 'x N' annotation from -m
}

defaults = {
	'input': 'AUTOBATCH', 
	'silent': False, 
	'opaque': False, 
	'multi_range': None,
	'color_exp': 'k', 
	'color_cal': 'r', 
	'color_pos': 'b', 
	'color_dif': 'g',
	'marker_exp_size': 6, 
	'marker_pos_size': 10, 
	'plot_size': (6, 4), 
	'dots_per_inch': 600, 
	'extension': 'svg', 
	'legend_text_exp': 'Observed', 
	'legend_text_cal': 'Calculated', 
	'legend_text_pos': 'Reflections', 
	'legend_text_dif': 'Difference', 
	'legend_columns': 1, 
	'size_axis_labels': 11, 
	'size_legend_labels': 12, 
	'size_tick_labels': 10, 
	'size_multiply_label': 12, 
	'x_label_text': r'$2\theta \quad / \quad°$', 
	'y_label_text': 'Intensity / a.u.', 
	'x_step_width': 5, 
	'vline_style': '-', 
	'vline_strength': 0.6
}

parser = argparse.ArgumentParser(description='Plots the result of a TOPAS Pawley fit using output files.')
parser.add_argument('-i', '--input', type=str, nargs='+', default=defaults['input'])
parser.add_argument('-s', '--silent', action='store_true', default=defaults['silent'])
parser.add_argument('-c', '--cell_info', action='store_true', help='Include unit cell parameter boxes on the plot.')
parser.add_argument('-m', '--multiply', nargs='+', help='Format: a,b,N. Use ,b,N or a,,N for limits.')
parser.add_argument('-x', '--extension', type=str, default=settings['extension'],
                    help='Output image format used with -s, e.g. svg, png, pdf (default: %(default)s).')
parser.add_argument('--qall', action='store_true',
                    help='Show all three fit-quality factors (R_wp, R_exp, chi) instead of R_wp alone.')
args = parser.parse_known_args()[0]
settings['extension'] = args.extension.lstrip('.').lower()

sgs_HM = {
	'1':   r'$P1$',
	'2':   r'$P\bar{1}$',
	'10':  r'$P2/m$',
	'11':  r'$P2_1/m$',
	'14':  r'$P2_1/c$',
	'19':  r'$P2_12_12_1$',
	'61':  r'$Pbca$',
	'110': r'$I4_1cd$',
	'148': r'$R\bar{3}$',
	'178': r'$P6_122$',
	'194': r'$P6_3/mmc$',
	'217': r'$I\bar{4}3d$',
}

# Aliases mapping raw HM strings (as they may appear in .out / .inp files,
# in any case and with/without underscore subscripts) to canonical sg-number keys.
# Lookup is done case-insensitively on a normalised form.
_hm_aliases = {
	'p1':         '1',
	'p-1':        '2',  'pbar1': '2',
	'p2/m':       '10',
	'p2_1/m':     '11', 'p21/m': '11',
	'p2_1/c':     '14', 'p21/c': '14',
	'p2_12_12_1': '19', 'p212121': '19',
	'pbca':       '61',
	'i4_1cd':     '110', 'i41cd': '110',
	'r-3':        '148', 'rbar3': '148', 'r3bar': '148',
	'p6_122':     '178', 'p6122': '178',
	'p6_3/mmc':   '194', 'p63/mmc': '194',
	'i-43d':      '217', 'ibar43d': '217', 'i4bar3d': '217',
}

def _norm_hm(s):
	"""Normalise an HM string for alias lookup: lowercase, strip whitespace."""
	return str(s).strip().lower() if s is not None else ''

def resolve_sg(raw):
	"""Given a raw space-group token (number string or HM symbol, with or without quotes),
	return (canonical_sg_number, latex_label).
	Either may be None / a best-effort fallback if the token isn't recognised."""
	if raw is None:
		return None, None
	s = str(raw).strip().strip('"').strip()
	if not s:
		return None, None
	if s.isdigit():
		return s, sgs_HM.get(s)
	key = _hm_aliases.get(_norm_hm(s))
	if key:
		return key, sgs_HM.get(key)
	# Unknown HM string — keep the raw token wrapped in LaTeX math mode as a fallback label.
	return None, f'${s}$'

# ==========================================
# COLOR TRANSFORM ENGINE
# ==========================================

def to_pastel(color, weight=settings['pastel_weight']):
	"""Blends an input color map key safely with pure white to generate a pastel background."""
	rgb = mcolors.to_rgb(color)
	return [(1.0 - weight) * c + weight for c in rgb]

# ==========================================
# CRYSTALLOGRAPHIC ROUNDING FUNCTION
# ==========================================

def cryst_round(mean_err):
	getcontext().prec = 32

	if mean_err is None:
		return None

	if '`_' in mean_err:
		mean, error = mean_err.split('`_', 1)
	elif '_' in mean_err:
		mean, error = mean_err.split('_', 1)
	else:
		return mean_err

	# TOPAS appends annotations after the numeric uncertainty (e.g. `_LIMIT_MAX_<v>`,
	# `_LIMIT_MIN_<v>`, `_SVD_ERR`, and other diagnostic suffixes). Since both the
	# mean and the uncertainty are always plain numbers, keep only the leading
	# numeric prefix of each side and discard whatever trails behind.
	_num = r'^\s*[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?'
	mn = re.match(_num, mean)
	en = re.match(_num, error)
	if mn:
		mean = mn.group(0)
	if en:
		error = en.group(0)

	# Salvage on partial garbage: if only the uncertainty side is malformed,
	# return the mean as a plain string (no `(esd)` annotation). If even the
	# mean can't parse, give up cleanly so the caller skips this row.
	try:
		mean  = Decimal(mean)
		error = Decimal(error)
	except (InvalidOperation, ValueError):
		try:
			return str(Decimal(mean))
		except (InvalidOperation, ValueError):
			return None
	
	# "Rule of 19": the bracketed esd is an integer from 2 to 19 — i.e. one
	# significant digit, or two when the leading digit is 1. An esd that would
	# round to 20 or more is shortened by one digit, and the mean is rounded to
	# the same decimal place: 15.4840(20) -> 15.484(2).
	error = abs(error)
	if error == 0:
		return format(mean, 'f')

	ndec = 1 - error.adjusted()  # decimal places that give the esd two significant digits
	while True:
		bracket = int(error.scaleb(ndec).to_integral_value(rounding=ROUND_HALF_UP))
		if bracket > 19:
			ndec -= 1
		elif bracket < 2:
			ndec += 1
		else:
			break

	mean_round = mean.quantize(Decimal(1).scaleb(-ndec), rounding=ROUND_HALF_UP)

	# Always render fixed-point — lattice parameters and cell volumes must never
	# appear in scientific notation (Decimal's str() flips to it for large
	# exponents). When the esd's significant digit sits left of the decimal
	# point (ndec < 0), pad it back out so the bracket stays unambiguous:
	# 6.6E+3 ± 7E+2 renders as 6600(700), not 6600(7).
	if ndec < 0:
		bracket *= 10 ** -ndec
	return '%s(%s)' % (format(mean_round, 'f'), bracket)

# ==========================================
# FILE WRANGLING FUNCTIONS
# ==========================================

def get_file_dicts():
	"""Discover TOPAS data file groups by recognising the data-type suffix at the end
	of each filename. Works regardless of whether the group uses `_pawley_NN_`,
	`_fit_NN_`, or no ident at all — the prefix up to the data-type marker becomes
	the group key."""
	all_group_dicts = {}
	# Recognised data-type suffixes. 2Th_Ip may optionally have an `_<sg>` suffix.
	data_idents = [
		settings['X_Yobs_ident'],
		settings['Out_X_Ycalc_ident'],
		settings['X_Difference_ident'],
		settings['2Th_Ip_ident'],
	]

	# Accept either a single extension string or any iterable of them
	exts = settings['pawley_ext']
	if isinstance(exts, str):
		exts = (exts,)
	all_data_files = [p for ext in exts for p in glob('*' + ext)]

	for f in all_data_files:
		base = os.path.splitext(os.path.basename(f))[0]
		for ident in data_idents:
			# Match: <prefix>_<ident>[_<extra>]   where <extra> is digits or HM-ish chars
			pat = rf'^(.*?)_({re.escape(ident)})(_[\w/+\-]+)?$'
			m = re.match(pat, base)
			if not m:
				continue
			prefix = m.group(1)
			ident_part = m.group(2) + (m.group(3) or '')
			all_group_dicts.setdefault(prefix, {})[ident_part] = f
			break
	return all_group_dicts


def sort_filegroup(group_dict):
	prios = {
		settings['X_Yobs_ident']: 0,
		settings['Out_X_Ycalc_ident']: 1,
		settings['2Th_Ip_ident']: 2,
		settings['X_Difference_ident']: 3,
	}
	
	def get_prio(tup):
		ident = tup[0]
		matched_keys = [key for key in prios if key in ident]
		return prios[matched_keys[0]] if matched_keys else 99

	return sorted(group_dict.items(), key=get_prio)


def get_x_y(path):
	filestring = Path(path).read_text()
	lines = [line.strip().split() for line in filestring.strip().split('\n') if line.strip()]
	x, y = np.array(lines, dtype=float).transpose()
	return x, y


def find_outfile_for_group(group_key, all_out_files=None):
	"""Locate the .out (TOPAS log) file that belongs to a data group.

	The .out filename is not always `{group_key}.out` — common variants:
	  - `{group_key}.out`                      (same ident as data files)
	  - `{base}.out`                            (no ident on .out)
	  - `{base}_pawley_NN.out` / `{base}_fit_NN.out`  (different ident than data)
	where `{base}` is `{group_key}` with any trailing `_pawley_NN` / `_fit_NN` stripped.

	Returns the matching filename or None if nothing is found.
	"""
	if all_out_files is None:
		all_out_files = glob('*.out')
	all_out_files = set(all_out_files)

	# 1. Exact ident match
	direct = f'{group_key}.out'
	if direct in all_out_files:
		return direct

	# 2. Strip trailing _pawley_NN / _fit_NN to get the bare "base" name
	base = re.sub(r'_(?:pawley|fit)_\d+$', '', group_key)
	if base != group_key and f'{base}.out' in all_out_files:
		return f'{base}.out'

	# 3. Any .out whose basename starts with base_<known-ident>
	for suffix in ('_pawley_01.out', '_pawley_02.out', '_fit_01.out', '_fit_02.out'):
		cand = f'{base}{suffix}'
		if cand in all_out_files:
			return cand

	# 4. Final fallback: any .out starting with `{base}_` (deterministic ordering)
	prefix_matches = sorted(f for f in all_out_files if f.startswith(f'{base}_'))
	if prefix_matches:
		return prefix_matches[0]

	return None

# ==========================================
# PARSING & DATA EXTRACTION
# ==========================================

def get_sgs_from_outfile(path):
	if not os.path.exists(path):
		return {}
		
	sgs = {}
	try:
		with open(path, 'r', encoding='utf-8') as inf:
			for line in inf:
				if 'Selected phases:' in line:
					matches = re.findall(r'"([^"]+)"\s*\((\d+)\)', line)
					for substance_name, sg_num in matches:
						sgs[sg_num] = substance_name
	except Exception as e:
		print(f"Warning: Error encountered reading metadata from {path}: {e}")
	return sgs


def get_outfile_info(path):
	"""Parse the fit-quality factors TOPAS writes on one line of the .out:
	`r_wp <v> r_exp <v> r_p <v> ... gof <v>`. Returns whichever of r_wp / r_exp /
	gof were found (keys are simply absent when the .out is missing or a factor
	can't be located), so the caller can skip annotations rather than invent a
	placeholder value. The `\\s+` after each key keeps the `_dash` variants
	(`r_wp_dash`, `r_exp_dash`) on the same line from being mistaken for the
	plain factors."""
	info = {}
	if not os.path.exists(path):
		return info

	try:
		filestring = Path(path).read_text()
		_float = r'([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)'
		for key in ('r_wp', 'r_exp', 'gof'):
			m = re.search(r'\b' + key + r'\s+' + _float, filestring)
			if m:
				info[key] = float(m.group(1))
	except Exception:
		pass
	return info


def _strip_topas_comments(text):
	"""Blank out TOPAS comment lines (lines whose first non-whitespace char is `'`).
	Returns the cleaned text with the same line count, so per-line regex matches
	preserve their semantic locations."""
	return '\n'.join(
		'' if re.match(r"\s*'", line) else line
		for line in text.split('\n')
	)


def _extract_cell_parts(label, raw_string):
	"""Round a raw `value` or `value`_uncertainty` token to crystallographic form and
	split into (label, main, bracket). Returns None if the token can't be parsed."""
	rounded = cryst_round(raw_string)
	if rounded is None:
		return None
	m = re.match(r"^([^(]+)(\(.*\))$", rounded)
	if m:
		return label, m.group(1), m.group(2)
	return label, rounded, ""


def get_unit_cell_info(path):
	"""Parse cell parameters, volume, and space group from a TOPAS .out file.
	Handles:
	  - Commented-out template lines (leading `'`) — ignored
	  - Single-line crystal-system macros: `Cubic(a)`, `Tetragonal(a, c)`,
	    `Orthorhombic(a, b, c)`, `Monoclinic(a, b, c, beta)`, `Triclinic(a, b, c, al, be, ga)`,
	    plus Hexagonal / Trigonal (same shape as Tetragonal)
	  - Loose per-line declarations: `a @ ...`, `b lpa ...`, `c @ ...`,
	    `al @ ...`, `be @ ...`, `ga @ ...` (also without the @)
	  - Volume keyword variants: `cell_volume` and `volume`
	  - Space-group syntaxes: quoted numeric `"61"`, quoted HM `"P63/mmc"`,
	    quoted lowercase `"p1"`, unquoted `Pbca`, hyphenated `R-3`

	Returns a list of (raw_sg_token, cell_data) tuples, one per phase, in the
	order the .out file declares them. raw_sg_token may be either a number string
	or an HM symbol; downstream callers should pass it through `resolve_sg`.
	"""
	if not os.path.exists(path):
		return []

	try:
		raw = Path(path).read_text(encoding='utf-8', errors='replace')
	except Exception as e:
		print(f"Warning: cannot read {path}: {e}")
		return []

	content = _strip_topas_comments(raw)

	# Each phase begins at `hkl_Is`. The preamble before the first marker is dropped.
	blocks = re.split(r'\bhkl_Is\b', content)
	phase_blocks = blocks[1:] if len(blocks) > 1 else [content]

	# Volumes (cell_volume or volume) collected in document order, assigned to phases
	# by index. They can appear inside or outside the hkl_Is block.
	all_vols = re.findall(r'\b(?:cell_volume|volume)\s+(\S+)', content)

	# Crystal-system macros recognised on a single line, e.g.  Cubic(@ 17.09...)
	macro_pat = re.compile(
		r'\b(Cubic|Tetragonal|Hexagonal|Trigonal|Orthorhombic|Monoclinic|Triclinic)'
		r'\s*\(([^)]+)\)'
	)

	# Loose per-line declarations. The leading word boundary keeps `a` from
	# matching `al`, `axial`, etc.; the (?:[a-zA-Z@]+\s+)* skips tokens like
	# `@` or `lpa` between the key and the numeric value.
	loose_param_keys = (('a', 'a'), ('b', 'b'), ('c', 'c'),
	                    ('alpha', 'al'), ('beta', 'be'), ('gamma', 'ga'))

	phases_data = []

	# Safety cap to avoid pathological inputs creating dozens of phantom phases
	for phase_idx, block in enumerate(phase_blocks[:5]):
		if not block.strip():
			continue

		# Space group: quoted form preferred (it can contain '/' and '-'),
		# otherwise the first whitespace-delimited token after `space_group`.
		sg_m = re.search(r'space_group\s+"([^"]+)"', block) \
		       or re.search(r'space_group\s+(\S+)', block)
		sg_raw = sg_m.group(1).strip() if sg_m else None

		cell_data = []

		# 1. Try the crystal-system macro form
		sys_m = macro_pat.search(block)
		if sys_m:
			param_str = sys_m.group(2)
			raw_tokens = [p.replace('@', '').strip() for p in param_str.split(',') if p.strip()]
			n_params = len(raw_tokens)
			if n_params == 1:
				labels = ['a']
			elif n_params == 2:
				labels = ['a', 'c']
			elif n_params == 3:
				labels = ['a', 'b', 'c']
			elif n_params == 4:
				labels = ['a', 'b', 'c', r'\beta']
			elif n_params == 6:
				labels = ['a', 'b', 'c', r'\alpha', r'\beta', r'\gamma']
			else:
				labels = [f'p{i+1}' for i in range(n_params)]
			for label, token in zip(labels, raw_tokens):
				try:
					parts = _extract_cell_parts(label, token)
				except Exception as e:
					print(f"Warning: skipping cell param '{label}' (couldn't parse '{token}'): {e}")
					parts = None
				if parts:
					cell_data.append(parts)
		else:
			# 2. Loose per-line declarations
			for label, key in loose_param_keys:
				pat = rf'^\s*{key}\b\s+(?:[a-zA-Z@]+\s+)*(\S+)'
				m = re.search(pat, block, re.MULTILINE)
				if m:
					try:
						parts = _extract_cell_parts(label, m.group(1))
					except Exception as e:
						print(f"Warning: skipping cell param '{label}' (couldn't parse '{m.group(1)}'): {e}")
						parts = None
					if parts:
						cell_data.append(parts)

		# 3. Volume (one per phase, by document order)
		if phase_idx < len(all_vols):
			try:
				parts = _extract_cell_parts('V', all_vols[phase_idx])
			except Exception as e:
				print(f"Warning: skipping volume (couldn't parse '{all_vols[phase_idx]}'): {e}")
				parts = None
			if parts:
				cell_data.append(parts)

		phases_data.append((sg_raw, cell_data))

	return phases_data

# ==========================================
# MATPLOTLIB COMPOSITION ARRAYS
# ==========================================

def stack_artists_vertically(ax, N_lines,
                             bragg_spacing=settings['bragg_spacing'],
                             top_clearance=settings['tick_top_clearance'],
                             bottom_clearance=settings['tick_bottom_clearance'],
                             y_tol_top=settings['y_tol_top'],
                             y_tol_bottom=settings['y_tol_bottom'],
                             d_difference=settings['difference_band_height'],
                             min_data_fraction=settings['min_data_fraction'],
                             diff_shift_factor=settings['diff_shift_factor']):
	"""
	Positions Bragg tick rows and the difference curve using axes-coordinate fractions
	so the layout is invariant to data amplitude (e.g. range multiplications).

	Layout from top to bottom (axes coords, 0 = bottom edge, 1 = top edge):
	  y_tol_top | data region | top_clearance | tick rows | bottom_clearance | difference | y_tol_bottom
	"""
	try:
		fig = ax.figure
		fig.canvas.draw()

		N_pos = N_lines - 3  # number of Bragg tick row artists
		if N_pos <= 0 or len(ax.lines) < N_lines:
			return

		# Physical marker height expressed as an axes-coordinate fraction
		bbox = ax.get_window_extent()
		if bbox.height == 0:
			return
		marker_height_axes = (settings['2Th_Ip_markersize'] / 72.0 * fig.dpi) / bbox.height

		# All N_pos tick rows, each one marker tall, separated by bragg_spacing
		total_ticks_frac = N_pos * marker_height_axes + max(0, N_pos - 1) * bragg_spacing

		# Actual data range from observed + calculated profiles
		y0 = ax.lines[0].get_ydata()
		y1 = ax.lines[1].get_ydata()
		y_data_min = min(y0.min(), y1.min())
		y_data_max = max(y0.max(), y1.max())
		dy_data = y_data_max - y_data_min
		if dy_data == 0:
			return

		# Difference curve amplitude (full peak-to-trough span). Used both to anchor
		# the curve and to grow d_difference if the natural amplitude wouldn't fit.
		# We read it before laying out so the band can be expanded in a single pass.
		diff_line = ax.lines[N_lines - 1]
		y_diff = diff_line.get_ydata()
		diff_amplitude = float(y_diff.max() - y_diff.min())

		# Everything except the data region and the diff band is fixed
		fixed_non_diff = (y_tol_top + y_tol_bottom + bottom_clearance
		                  + total_ticks_frac + top_clearance)

		# Auto-expand d_difference so the natural diff amplitude fits inside
		# `diff_shift_factor` of the band (e.g. 0.8 leaves a 20% safety margin).
		# Derivation: we require D * dy_data / (1 - fixed_non_diff - D) ≥ diff_amplitude / sf
		# Solving for the minimum D yields the expression below.
		if diff_amplitude > 0 and 0 < diff_shift_factor <= 1:
			sf = diff_shift_factor
			denom = dy_data * sf + diff_amplitude
			if denom > 0:
				d_diff_needed = diff_amplitude * (1.0 - fixed_non_diff) / denom
				d_difference = max(d_difference, d_diff_needed)

		# Data region: whatever's left after reserving every other band
		d_calc_exp = 1.0 - (fixed_non_diff + d_difference)
		d_calc_exp = max(d_calc_exp, min_data_fraction)

		# Compute ylim so data fills exactly d_calc_exp of the axes height,
		# with y_data_max sitting at axes coord (1 - y_tol_top).
		ylim_range = dy_data / d_calc_exp
		ylim_max   = y_data_max + y_tol_top * ylim_range
		ylim_min   = ylim_max - ylim_range

		# Axes-coord → data-coord converter
		def ac2dc(ac):
			return ylim_min + ac * ylim_range

		# Axes coord where the data floor (y_data_min) sits
		data_bottom_axes = 1.0 - y_tol_top - d_calc_exp

		# Place each Bragg tick row below the data floor
		first_tick_center_axes = data_bottom_axes - top_clearance - marker_height_axes / 2
		for i in range(N_pos):
			tick_center_axes = first_tick_center_axes - i * (marker_height_axes + bragg_spacing)
			tick_y_data = ac2dc(tick_center_axes)
			line_idx = i + 2
			ax.lines[line_idx].set_ydata(
				np.full_like(ax.lines[line_idx].get_xdata(), tick_y_data)
			)

		# Anchor the difference curve at the centre of d_difference.
		# Using the median (rather than the geometric midpoint) is robust against
		# the occasional large spike. Note: no multiplier here — the previous
		# `* diff_shift_factor` was an operator-precedence bug that shifted the
		# curve upward (toward the Bragg ticks) by ~20% of |ac2dc(baseline)|.
		# Use N_lines - 1 (not -1) because process_multiplication may have appended
		# axvline objects to ax.lines after the original data lines were plotted.
		diff_baseline_axes = y_tol_bottom + d_difference / 2
		shift = float(np.median(y_diff)) - ac2dc(diff_baseline_axes)
		diff_line.set_ydata(y_diff - shift)

		ax.set_ylim(ylim_min, ylim_max)

	except Exception as e:
		print(f"Warning: Skipping vertical scaling alignment due to formatting error: {e}")


# ==========================================
# MULTIPLICATION LOGIC
# ==========================================

def process_multiplication(ax, multi_args):
	# exp (line 0), calc (line 1), diff (last line)
	target_indices = [0, 1, -1]
	
	for arg in multi_args:
		try:
			parts = arg.split(',')
			a_str, b_str, n_str = parts[0], parts[1], parts[2]
			
			# Determine x-range
			x_min_data = ax.lines[1].get_xdata().min()
			x_max_data = ax.lines[1].get_xdata().max()
			
			a = float(a_str) if a_str else x_min_data
			b = float(b_str) if b_str else x_max_data
			n = float(n_str)
			
			# Apply to exp, calc, and diff
			for idx in target_indices:
				line = ax.lines[idx]
				x, y = line.get_data()
				mask = (x >= a) & (x <= b)
				y[mask] *= n
				line.set_data(x, y)
			
			for x_val in [a, b]:
				ax.axvline(x=x_val, color='k', linestyle='--', linewidth=0.8, alpha=0.5,
				           label='_nolegend_')
			
			label_text = f"x {n}"
			# Dynamic placement: anchor the label to the boundary of [a, b] closest to
			# the plot centre and align horizontally inward, so the text sits inside
			# the multiplied range and stays clear of the upper-right legend.
			# Clip to the data range so user-supplied bounds outside the data don't
			# push the label off the axes.
			a_clip = max(a, x_min_data)
			b_clip = min(b, x_max_data)
			plot_mid = (x_min_data + x_max_data) / 2
			range_mid = (a_clip + b_clip) / 2
			dead_zone = (x_max_data - x_min_data) * 0.1

			if range_mid < plot_mid - dead_zone:
				text_x, ha = b_clip, 'right'
			elif range_mid > plot_mid + dead_zone:
				text_x, ha = a_clip, 'left'
			else:
				text_x, ha = range_mid, 'center'

			# If the anchor is too close to a plot edge for the ha direction to fit
			# the text inside the axes, flip the alignment.
			anchor_ax = (text_x - x_min_data) / (x_max_data - x_min_data)
			if ha == 'left' and anchor_ax > 0.85:
				ha = 'right'
			elif ha == 'right' and anchor_ax < 0.15:
				ha = 'left'

			trans = blended_transform_factory(ax.transData, ax.transAxes)
			ax.text(text_x, settings['multiply_label_y'], label_text,
			        ha=ha, va='top', fontsize=10, transform=trans)
			
		except (ValueError, IndexError) as e:
			print(f"Skipping invalid multiplication argument {arg}: {e}")


def style(ax):
	size_axis_labels = getattr(args, 'size_axis_labels', defaults['size_axis_labels'])
	size_tick_labels = getattr(args, 'size_tick_labels', defaults['size_tick_labels'])

	ax.set_xlabel(r'$2\theta \quad / \quad \mathrm{^\circ}$', size=size_axis_labels, labelpad=7)
	ax.set_ylabel(r'$\mathrm{Intensity} \quad / \quad \mathrm{a.u.}$', size=size_axis_labels, labelpad=7)
	ax.set_yticks([])
	ax.xaxis.set_minor_locator(AutoMinorLocator())
	ax.tick_params(axis='both', which='both', labelsize=size_tick_labels, direction='in')
	
	if len(ax.lines) > 1:
		xmin = ax.lines[1].get_xdata().min()
		xmax = ax.lines[1].get_xdata().max()
		ax.set_xlim(xmin, xmax)


def add_legend(ax):
	visible = [l for l in ax.lines if not l.get_label().startswith('_')]
	leg_lines = [Line2D([0], [0], marker=l.get_marker(), ms=l.get_ms(), ls=l.get_ls(), c=l.get_color(), lw=1) for l in visible]
	leg_labs = [l.get_label() for l in visible]
	ax.legend(leg_lines, leg_labs, fontsize=settings['legend_fontsize'], frameon=True, loc='upper right')


def add_quality(ax, info, show_all=False):
	"""Annotate the fit-quality factor(s) in the bottom-right corner.

	Default: the weighted profile R-factor R_wp alone (falling back to chi if
	R_wp wasn't parsed). With ``show_all`` (the --qall flag): R_wp, R_exp and chi
	on one line, in the same italic style. R-factors carry a `%` since TOPAS
	reports them as percentages; chi is dimensionless. Any factor that wasn't
	found in the .out is silently dropped, and if none are available no text is
	drawn at all."""
	if show_all:
		wanted = [('R_{wp}', 'r_wp', r'\%'), ('R_{exp}', 'r_exp', r'\%'), (r'\chi', 'gof', '')]
	elif info.get('r_wp') is not None:
		wanted = [('R_{wp}', 'r_wp', r'\%')]
	else:
		wanted = [(r'\chi', 'gof', '')]

	parts = []
	for symbol, key, unit in wanted:
		val = info.get(key)
		if val is not None:
			parts.append(r'%s = %.2f%s' % (symbol, val, unit))
	if not parts:
		return

	text = '$' + r',\ '.join(parts) + '$'
	ax.text(.99, .01, text, ha='right', va='bottom', size=12, style='italic', transform=ax.transAxes)


def _renderer(fig):
	"""Best-effort renderer fetch that works across the Agg-family backends."""
	try:
		return fig.canvas.get_renderer()
	except AttributeError:
		fig.canvas.draw()
		return fig.canvas.get_renderer()


def _frac_bbox(artist, ax, renderer):
	"""Window extent of an artist (its rounded bbox patch, if any) in axes-fraction coords."""
	getter = getattr(artist, 'get_bbox_patch', None)
	target = getter() if getter and getter() is not None else artist
	bb = target.get_window_extent(renderer)
	(x0, y0), (x1, y1) = ax.transAxes.inverted().transform([(bb.x0, bb.y0), (bb.x1, bb.y1)])
	return x0, y0, x1, y1


def _data_envelope_frac(ax):
	"""The observed + calculated traces expressed in axes-fraction coordinates.
	transLimits maps data → axes fraction purely from the view limits, so the result
	is independent of figure size and can be computed once."""
	chunks = []
	for ln in ax.lines[:2]:  # plot order: [0] observed, [1] calculated
		xd = np.asarray(ln.get_xdata(), dtype=float)
		yd = np.asarray(ln.get_ydata(), dtype=float)
		if xd.size:
			chunks.append(ax.transLimits.transform(np.column_stack([xd, yd])))
	return np.vstack(chunks) if chunks else np.empty((0, 2))


def _envelope_ceiling(env, x_left, x_right):
	"""Highest data y-fraction within the axes-fraction x-span [x_left, x_right]."""
	if env.size == 0:
		return 0.0
	mask = (env[:, 0] >= x_left) & (env[:, 0] <= x_right)
	return float(env[mask, 1].max()) if mask.any() else 0.0


def _measure_boxes(boxes, fontsize, ax, fig):
	"""Set all boxes to `fontsize`, draw once, and return per-box
	(width, height, off_right, off_top) in axes fraction. The offsets are the gap
	between the text anchor and the patch's right/top edges, so a caller can place
	the patch's corner exactly by anchoring at (target - offset). All four values
	depend only on the font size and the figure's physical size, not on position."""
	for b in boxes:
		b.set_fontsize(fontsize)
	fig.canvas.draw()
	renderer = _renderer(fig)
	dims = []
	for b in boxes:
		x0, y0, x1, y1 = _frac_bbox(b, ax, renderer)
		ax_x, ax_y = b.get_position()
		dims.append((x1 - x0, y1 - y0, x1 - ax_x, y1 - ax_y))
	return dims


def _apply_vertical(boxes, dims, top, right, vgap):
	"""Right-aligned single column, stacked downward from `top`."""
	y = top
	for b, (w, h, off_r, off_t) in zip(boxes, dims):
		b.set_position((right - off_r, y - off_t))
		y -= h + vgap


def _apply_horizontal(boxes, dims, top, right, hgap):
	"""Right-aligned single row, top edges level with `top`."""
	r = right
	for b, (w, h, off_r, off_t) in zip(boxes, dims):
		b.set_position((r - off_r, top - off_t))
		r -= w + hgap


def add_unit_cell_boxes(ax, phases_cell_info, phase_colors, x=settings['box_x']):
	"""Place one info box per phase, coloured to match its Bragg reflections, packed
	tightly under the legend. The arrangement adapts to the space actually available
	between the legend and the plotted data, following this hierarchy:

	  1. vertical stack (one column) at the default font size;
	  2. else a horizontal row (side by side) at the default font size;
	  3. else shrink the font (retrying vertical, then horizontal) down to
	     `box_fontsize_min`;
	  4. else accept overlap, using whichever candidate descends least into the data.

	The layout is measured against the live renderer and recomputed on every resize,
	so the interactive window and the saved file agree regardless of window size."""
	if not phases_cell_info:
		return
	fig = ax.figure

	# Build the text artists once; positions and font size are set later by relayout.
	boxes = []
	for idx, (sg_num, cell_info) in enumerate(phases_cell_info):
		if not cell_info:
			continue
		color = phase_colors[idx] if idx < len(phase_colors) \
			else settings['2Th_Ip_colors'][idx % len(settings['2Th_Ip_colors'])]
		max_val_len = max(len(mv + brk) for _, mv, brk in cell_info)
		lines = [f"${lbl}$: {(mv + brk):>{max_val_len}}" for lbl, mv, brk in cell_info]
		boxes.append(ax.text(
			settings['box_x'], settings['box_y_top'], "\n".join(lines),
			transform=ax.transAxes, family='monospace',
			fontsize=settings['box_fontsize'], ha='right', va='top',
			multialignment='left',
			bbox=dict(boxstyle=f"round,pad={settings['box_pad']}",
					  facecolor=to_pastel(color), edgecolor='none')))
	if not boxes:
		return

	env = _data_envelope_frac(ax)
	clearance = settings['box_data_clearance']
	vgap = settings['box_gap']
	hgap = settings['box_col_gap']

	def relayout():
		renderer = _renderer(fig)

		# Anchor right below the legend; fall back to fixed coords if it's absent.
		leg = ax.get_legend()
		if leg is not None:
			lx0, ly0, lx1, ly1 = _frac_bbox(leg, ax, renderer)
			top, right = ly0 - settings['box_legend_gap'], lx1
		else:
			top, right = settings['box_y_top'], settings['box_x']

		best = None  # (arrangement, fontsize, block_bottom) — least-intrusive fallback
		for fs in range(settings['box_fontsize'], settings['box_fontsize_min'] - 1, -1):
			dims = _measure_boxes(boxes, fs, ax, fig)
			ws = [d[0] for d in dims]
			hs = [d[1] for d in dims]

			# Vertical: one right-aligned column.
			v_left = right - max(ws)
			v_bottom = top - (sum(hs) + vgap * (len(boxes) - 1))
			v_fits = v_left >= 0 and v_bottom >= _envelope_ceiling(env, v_left, right) + clearance

			# Horizontal: one right-aligned row.
			h_left = right - (sum(ws) + hgap * (len(boxes) - 1))
			h_bottom = top - max(hs)
			h_fits = h_left >= 0 and h_bottom >= _envelope_ceiling(env, h_left, right) + clearance

			if v_fits:
				_apply_vertical(boxes, dims, top, right, vgap)
				return
			if h_fits:
				_apply_horizontal(boxes, dims, top, right, hgap)
				return

			for arrange, bottom in (('V', v_bottom), ('H', h_bottom)):
				if best is None or bottom > best[2]:
					best = (arrange, fs, bottom)

		# Nothing fit at any size — accept overlap with the least-descending candidate.
		arrange, fs, _ = best
		dims = _measure_boxes(boxes, fs, ax, fig)
		(_apply_vertical if arrange == 'V' else _apply_horizontal)(
			boxes, dims, top, right, vgap if arrange == 'V' else hgap)

	relayout()

	# Re-pack when the interactive window is resized so the look stays consistent.
	state = {'busy': False}

	def _on_resize(event):
		if state['busy']:
			return
		state['busy'] = True
		try:
			relayout()
		finally:
			state['busy'] = False
		fig.canvas.draw_idle()

	fig.canvas.mpl_connect('resize_event', _on_resize)

# ==========================================
# MAIN ROUTINE EXECUTION
# ==========================================

def _pick_by_ident(sorted_filegroup, ident_substr):
	"""First (ident, path) in the group whose ident contains the substring, or None."""
	for ident, path in sorted_filegroup:
		if ident_substr in ident:
			return ident, path
	return None


def _tick_sg_token(ident):
	"""Return the suffix after `2Th_Ip_` in a tick-file ident, or None if absent."""
	m = re.search(rf'^{re.escape(settings["2Th_Ip_ident"])}_(\S+)$', ident)
	return m.group(1).strip() if m else None


def _legend_label_for(canon_sg, hm_label, substance):
	"""Build a legend entry from whatever pieces are available, gracefully degrading
	when the substance is unknown."""
	fmt = settings['use_sg_format'] if settings['use_sg_from_outfile'] else None
	if fmt == 'SUBSTANCE' and substance:
		return substance
	if fmt == 'NUMBER' and canon_sg:
		return canon_sg
	if fmt == 'HERMANN-MAUGUIN' and hm_label:
		return hm_label
	if fmt == 'HERMANN-MAUGUIN+SUBSTANCE':
		if hm_label and substance:
			return f"{hm_label} ({substance})"
		if hm_label:
			return hm_label
		if substance:
			return substance
	# Last-ditch fallbacks
	if hm_label:
		return hm_label
	if canon_sg:
		return canon_sg
	return settings['2Th_Ip_label'] if isinstance(settings['2Th_Ip_label'], str) else 'Reflections'


file_dicts = get_file_dicts()
all_out_files = glob('*.out')

for group_name in file_dicts:
	group_dict = file_dicts[group_name]
	if len(group_dict) < 3:
		continue

	sorted_filegroup = sort_filegroup(group_dict)

	# Resolve files by type rather than positional index, so groups with missing
	# pieces are detected explicitly instead of silently mis-pairing.
	exp_file  = _pick_by_ident(sorted_filegroup, settings['X_Yobs_ident'])
	calc_file = _pick_by_ident(sorted_filegroup, settings['Out_X_Ycalc_ident'])
	dif_file  = _pick_by_ident(sorted_filegroup, settings['X_Difference_ident'])
	pos_files = [t for t in sorted_filegroup if settings['2Th_Ip_ident'] in t[0]]

	missing = [name for name, val in
	           (('X_Yobs', exp_file), ('Out_X_Ycalc', calc_file),
	            ('X_Difference', dif_file), ('2Th_Ip', pos_files)) if not val]
	if missing:
		print(f"Skipping group '{group_name}': missing file type(s) {missing}.")
		continue

	# Locate the .out — may be at a different basename than the data files,
	# or absent entirely (in which case metadata silently disappears but the
	# core obs/calc/ticks/difference still plot).
	outfile_path = find_outfile_for_group(group_name, all_out_files) or ''
	outfile_info     = get_outfile_info(outfile_path)
	phases_cell_info = get_unit_cell_info(outfile_path)   # [(sg_raw, cell_data), ...]
	substances       = get_sgs_from_outfile(outfile_path) # {sg_num: name} or {}

	# Pre-compute the canonical sg-number for each phase (for sg↔tick matching)
	phase_canon_sgs = [resolve_sg(sg_raw)[0] for sg_raw, _ in phases_cell_info]

	# Build per-tick metadata: filename sg (if any), positional fallback to phase
	# sg, then color and legend label.
	tick_meta = []
	for i, (ident, path) in enumerate(pos_files):
		sg_token = _tick_sg_token(ident)
		canon_sg, hm_label = resolve_sg(sg_token)
		# If the tick filename carried no sg, fall back to the i-th phase's sg
		# (best-effort: relies on declaration order matching).
		if not canon_sg and i < len(phase_canon_sgs) and phase_canon_sgs[i]:
			canon_sg = phase_canon_sgs[i]
			hm_label = sgs_HM.get(canon_sg, hm_label)
		substance = substances.get(canon_sg) if canon_sg else None
		tick_meta.append({
			'idx':    i,
			'path':   path,
			'color':  settings['2Th_Ip_colors'][i % len(settings['2Th_Ip_colors'])],
			'sg_num': canon_sg,
			'label':  _legend_label_for(canon_sg, hm_label, substance),
		})

	# Pair phases (from .out) to ticks in legend order. sg-match first, then take
	# the next remaining phase positionally. Boxes therefore appear in the same
	# top-to-bottom order as the legend's Bragg entries.
	phases_remaining = list(phases_cell_info)
	ordered_phases     = []
	ordered_box_colors = []
	for tick in tick_meta:
		matched, matched_j = None, None
		if tick['sg_num']:
			for j, (sg_raw, cell_data) in enumerate(phases_remaining):
				canon, _ = resolve_sg(sg_raw)
				if canon == tick['sg_num']:
					matched, matched_j = (sg_raw, cell_data), j
					break
		if matched is None and phases_remaining:
			matched, matched_j = phases_remaining[0], 0
		if matched is not None:
			ordered_phases.append(matched)
			ordered_box_colors.append(tick['color'])
			phases_remaining.pop(matched_j)

	# Tick legend labels in pos_files order
	settings['2Th_Ip_label'] = [tick['label'] for tick in tick_meta]

	# ---- Plot the core artists ----
	fig, ax = plt.subplots(figsize=settings['figsize'], layout='constrained')

	# Trace 1: Observed
	x, y = get_x_y(exp_file[1])
	ax.plot(x, y, ls='', marker='x', mew=.6, color=settings['X_Yobs_color'],
	        label=settings['X_Yobs_label'], ms=settings['X_Yobs_markersize'])

	# Trace 2: Calculated
	x, y = get_x_y(calc_file[1])
	ax.plot(x, y, ls='-', lw=1, marker='', color=settings['Out_X_Ycalc_color'],
	        label=settings['Out_X_Ycalc_label'], ms=settings['Out_X_Ycalc_markersize'])

	# Trace 3+: Bragg tick rows (one per phase)
	for i, (ident, path) in enumerate(pos_files):
		x, y = get_x_y(path)
		ax.plot(x, np.zeros_like(x), ls='', marker='|', mew=.6,
		        color=tick_meta[i]['color'], label=tick_meta[i]['label'],
		        ms=settings['2Th_Ip_markersize'])

	# Trace N: Difference
	x, y = get_x_y(dif_file[1])
	ax.plot(x, y, ls='-', marker='', lw=.3, color=settings['X_Difference_color'],
	        label=settings['X_Difference_label'], ms=settings['X_Difference_markersize'])

	# Layout, legend, optional decorations
	N_files = 2 + len(pos_files) + 1  # exp + calc + Bragg rows + diff
	if args.multiply:
		process_multiplication(ax, args.multiply)
	stack_artists_vertically(ax, N_files)
	add_legend(ax)

	if settings['show_info']:
		add_quality(ax, outfile_info, show_all=args.qall)
	if args.cell_info:
		add_unit_cell_boxes(ax, ordered_phases, ordered_box_colors)
	style(ax)

	if args.silent:
		outfile_name = f"{group_name}.{settings['extension']}"
		plt.savefig(outfile_name, dpi=settings['dpi'], bbox_inches='tight',
		            transparent=settings['transparent'])
		plt.close(fig)
	else:
		plt.show()