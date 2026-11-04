import matplotlib.pyplot as plt

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--brick", type=int, default=22)
parser.add_argument("--logdir", default="/afs/cern.ch/user/l/luarella/work-luarella/public/condor_logs/matching_sweep")
parser.add_argument("--xmin", type=float)
parser.add_argument("--xmax", type=float)
parser.add_argument("--ymin", type=float)
parser.add_argument("--ymax", type=float)
parser.add_argument("--xstep", type=float, default=0.01)
parser.add_argument("--ystep", type=float, default=0.01)
args = parser.parse_args()

# arrow positions, directions, errors, and amplitudes
l_x_pos = []
l_y_pos = []
l_x_dir = []
l_y_dir = []
l_x_err = []
l_y_err = []
l_x_amp = []
l_y_amp = []
l_x_bkg_amp = []
l_y_bkg_amp = []

xi = args.xmin
xf = args.xmax
yi = args.ymin
yf = args.ymax
xstep = args.xstep
ystep = args.ystep

x = xi
while x <= xf:
    y = yi
    while y <= yf:
        y = round(y,2)
        """The log file contains a bunch of stuff, but we are
        interested in the result of fits on the x and y axes. This
        loop tries to grab both fits. They come in sequence, first the
        fit in the x axis, then in y.

        Each of the fits consists of the sum of two gaussians, one
        centered around zero (due to random background), the other has
        the center as a free parameter (which indicates the
        "direction" of the arrow, which represents a bias).

        The parameters are
        - p0: amplitude of the bias gaussian
        - p1: center of the bias gaussian
        - p2: spread of the bias gaussian (not used in this macro)
        - p3: the amplitude of the background gaussian
        - p4: spread of the backgrond gaussian

        """
        filename = f"{args.logdir}/b{args.brick}_{x:.2f}_{y:.2f}.log"
        try:
            with open(filename) as f:
                print(f"{x} {y} opened")
                is_y=-1 # we use y=-1 as "neither x nor y", y=0 as x=true, and y=1 as y=true, x=false
                x_invalid = False
                y_invalid = False
                x_dir = 0.0
                y_dir = 0.0
                x_err = 0.0
                y_err = 0.0
                for line in f.readlines():
                    # line full of asterisks delimites minuit fit result
                    if "*************************************" in line:
                        is_y += 1
                    if "Invalid" in line:
                        is_y -= 1 # the invalid flag adds another line of asterisks
                        if is_y:
                            y_invalid = True
                            print(f"{x} {y}: Y invalid")
                        else:
                            x_invalid = True
                            print(f"{x} {y}: X invalid")
                    if is_y == -1:
                        continue
                    
                    split = line.split()
                    #print(f"{x} {y} split {split}")
                    if split[0] == "p0": # bias gaussian amplitude
                        if is_y:
                            y_amp = float(split[2])
                        else:
                            x_amp = float(split[2])
                    elif split[0] == "p1": # bias gaussian amplitude
                        if is_y:
                            y_dir = float(split[2])
                            y_err = float(split[4])
                        else:
                            x_dir = float(split[2])
                            x_err = float(split[4])
                    elif split[0] == "p3":
                        if is_y:
                            y_bkg_amp = float(split[2])
                            break # last thing to grab!
                        else:
                            x_bkg_amp = float(split[2])
                    else:
                        continue

                if x_err == 0.0:
                    print(f"WARNING: {x} {y}: x error {x_err}. setting to 0.0001")
                    x_err = 0.0001
                if y_err == 0.0:
                    print(f"WARNING: {x} {y}: y error {y_err}. setting to 0.0001")
                    y_err = 0.0001
                if not x_invalid and not y_invalid and abs(x_dir) < 0.5 and abs(y_dir) < 0.5:
                    l_x_pos.append(x)
                    l_x_dir.append(x_dir)
                    l_y_pos.append(y)
                    l_y_dir.append(y_dir)
                    l_x_amp.append(x_amp)
                    l_x_bkg_amp.append(x_bkg_amp)
                    l_y_amp.append(y_amp)
                    l_y_bkg_amp.append(y_bkg_amp)
                    l_x_err.append(x_err)
                    l_y_err.append(y_err)
        except:
            print(f"no file {filename}")
        y += ystep
    x += xstep
    x = round(x,2)

xlabel = f"Emulsion global translation in $x$ (cm)"
ylabel = f"Emulsion global translation in $y$ (cm)"

fig,ax = plt.subplots(figsize = (12,12))
ax.quiver(l_x_pos, l_y_pos, l_x_dir, l_y_dir)
ax.set_title("Extracted bias")
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
plt.savefig(f"img/quiver_dir_b{args.brick}.pdf")

fig,ax = plt.subplots(figsize = (12,12))
ax.quiver(l_x_pos, l_y_pos, l_x_err, l_y_err)
ax.set_title("Extracted bias errors")
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
plt.savefig(f"img/quiver_err_b{args.brick}.pdf")

# weighted per amplitudes
l_x_dir_w = [ x*(sig/bkg) for x, sig, bkg in zip(l_x_dir, l_x_amp, l_x_bkg_amp) ]
l_y_dir_w = [ y*(sig/bkg) for y, sig, bkg in zip(l_y_dir, l_y_amp, l_y_bkg_amp) ]

from math import sqrt
threshold = 2.0
filtered = []
for x, y, dx, dy in zip(l_x_pos, l_y_pos, l_x_dir_w, l_y_dir_w):
    mag = sqrt(dx*dx+dy*dy)
    if mag <= threshold:
        filtered.append((x,y,dx,dy))
if filtered:
    l_x_pos_f, l_y_pos_f, l_x_dir_w_f, l_y_dir_w_f = zip(*filtered)
else:
    l_x_pos_f, l_y_pos_f, l_x_dir_w_f, l_y_dir_w_f = [], [], [], []

fig,ax = plt.subplots(figsize = (12,12))
ax.quiver(l_x_pos_f, l_y_pos_f, l_x_dir_w_f, l_y_dir_w_f)
ax.set_title("Extracted bias weighted by amplitude of signal/background")
ax.set_xlabel(xlabel)
ax.set_ylabel(ylabel)
plt.savefig(f"img/quiver_dirw_b{args.brick}.pdf")
