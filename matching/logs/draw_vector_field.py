import matplotlib.pyplot as plt

from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument("--data", action='store_true')
parser.add_argument("--brick", type=int, default=22)
parser.add_argument("--logdir", default="/afs/cern.ch/user/l/luarella/work-luarella/public/condor_logs/matching_sweep")
parser.add_argument("--minx", type=float)
parser.add_argument("--maxx", type=float)
parser.add_argument("--miny", type=float)
parser.add_argument("--maxy", type=float)
parser.add_argument("--step", type=float, default=0.01)
args = parser.parse_args()

# arrow pos and dir
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

if args.data:
    #xi = 0.2
    #xf = 0.8
    xi = -0.10
    xf = 0.10
    xstep = 0.01

    #yi = -0.8
    #yf = -0.2
    yi = 0.05
    yf = 0.20
    ystep = 0.01
else:
    xi = -0.3
    xf = 0.3
    xstep = 0.01

    yi = -0.1
    yf =  0.1
    ystep = 0.01

x = xi
while x <= xf:
    y = yi
    while y <= yf:
        #print(x, y)
        y = round(y,2)
        if args.data:
            filename = f"data_s2only/{x:.2f}_{y:.2f}.log"
        else:
            filename = f"mc/{x:.2f}_{y:.2f}.log"
        #print(filename)
        try:
            #print(f"{x} {y} trying")
            with open(filename) as f:
                print(f"{x} {y} opened")
                is_y=-1
                x_invalid = False
                y_invalid = False
                x_dir = 0.0
                y_dir = 0.0
                x_err = 0.0
                y_err = 0.0
                for line in f.readlines():
                    if "*************************************" in line:
                        #print(f"{x} {y} found ast, isy {is_y}")
                        is_y += 1
                    if "Invalid" in line:
                        is_y -= 1 # the invalid flagg adds another line of asterisks
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
                    if split[0] == "p0":
                        if is_y:
                            y_amp = float(split[2])
                        else:
                            x_amp = float(split[2])
                    elif split[0] == "p1":
                        if is_y:
                            y_dir = float(split[2])
                            y_err = float(split[4])
                        else:
                            x_dir = float(split[2])
                            x_err = float(split[4])
                    elif split[0] == "p3":
                        if is_y:
                            y_bkg_amp = float(split[2])
                            break # last thing to grab!!!!!!!!!!!
                        else:
                            x_bkg_amp = float(split[2])
                    else:
                        continue
                #print(f"{x} {y} finished readlines")
                    
                if x == 0.47 and y == -0.73:
                    print("my event:")
                    print(x_dir, x_err, x_amp, x_bkg_amp)
                    print(y_dir, y_err, y_amp, y_bkg_amp)
                #print(f"{x} {y} append pos x")

                if not x_invalid and not y_invalid and abs(x_dir) < 0.5 and abs(y_dir) < 0.5:
                    l_x_pos.append(x)
                    l_x_dir.append(x_dir)
                    l_y_pos.append(y)
                    l_y_dir.append(y_dir)
                    l_x_amp.append(x_amp)
                    l_x_bkg_amp.append(x_bkg_amp)
                    l_y_amp.append(y_amp)
                    l_y_bkg_amp.append(y_bkg_amp)
                if x_err == 0.0:
                    #print(f"WARNING: {x} {y}: x error {x_err}. setting to 0.0001")
                    x_err = 0.0001
                if y_err == 0.0:
                    #print(f"WARNING: {x} {y}: y error {y_err}. setting to 0.0001")
                    y_err = 0.0001
                if not x_invalid and not y_invalid and abs(x_dir) < 0.5 and abs(y_dir) < 0.5:
                    #print(f"{x} {y} append x")
                    l_x_err.append(x_err)
                    #print(f"{x} {y} append y")
                    l_y_err.append(y_err)
                #print(f"{x} {y} finished")
        except:
            print(f"no file {filename}")
        y += ystep
    x += xstep
    x = round(x,2)

fig,ax = plt.subplots(figsize = (12,12))
ax.quiver(l_x_pos, l_y_pos, l_x_dir, l_y_dir)
ax.set_title("Averages")
if args.data:
    plt.savefig("img/quiver_dir_data_s2only.pdf")
else:
    plt.savefig("img/quiver_dir_mc.pdf")



fig,ax = plt.subplots(figsize = (12,12))
ax.quiver(l_x_pos, l_y_pos, l_x_err, l_y_err)
ax.set_title("Errors")
if args.data:
    plt.savefig("img/quiver_err_data_s2only.pdf")
else:
    plt.savefig("img/quiver_err_mc.pdf")


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
ax.set_title("weighted")
if args.data:
    plt.savefig("img/quiver_dirw_data_s2only.pdf")
else:
    plt.savefig("img/quiver_dirw_mc.pdf")

