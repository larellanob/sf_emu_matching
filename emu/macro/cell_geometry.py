import ROOT

def emu_to_sf(emux, emuy, brick):
    # translations (in cm
    x_t = -45.9
    y_t =  15.2
    if brick%10 == 1 or brick%10 == 3:
        x_t += 19.5
    if brick%10 == 3 or brick%10 == 4:
        y_t += 19.5
    return (emux/10000.)+x_t, (emuy/10000.)+y_t

def get_cell_borders(cellx,celly,brick):
    # all in emu coordinate system (um, 0 to 18000)
    width=10000 # cell width 1cm
    extra=0
    off_x, off_y = 0, 0
    '''
    off_x, off_y = 0, -500
    if brick%10 == 1 or brick%10 == 3:
        off_x = 195000
    if brick%10 == 3 or brick%10 == 4:
        off_y = 195000
    '''
    xi = off_x + cellx*width - width/2 - extra
    xf = off_x + cellx*width + width/2 + extra
    yi = off_y + celly*width - width/2 - extra
    yf = off_y + celly*width + width/2 + extra
    return xi, xf, yi, yf
