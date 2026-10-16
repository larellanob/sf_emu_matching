import ROOT

def get_brick_translation(brick, data):
    # returns translation in cm
    x_t = -45.9
    y_t =  15.2
    if data:
        x_t += 19.5
        y_t += 00.0
    return 

def emu_to_sf(emux, emuy, x_t, y_t):
    return (emux/10000.)+x_t, (emuy/10000.)+y_t

def get_cell_borders(x,y,brick,data):
    # all in emu coordinate system (um, 0 to 18000)
    width=10000
    extra=0
    off_x, off_y = 0, -500
    if brick%10 == 3 or brick%10 == 1:
        off_x = 195000
    if brick%10 == 3 or brick%10 == 4:
        off_y = 195000
    if data:
        off_x = 0
        off_y = 0
    xi = off_x + x*width - width/2 - extra
    xf = off_x + x*width + width/2 + extra
    yi = off_y + y*width - width/2 - extra
    yf = off_y + y*width + width/2 + extra
    return xi, xf, yi, yf

