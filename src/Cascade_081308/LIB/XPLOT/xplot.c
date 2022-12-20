/*

  graphic interface XPLOT to perform basic graphic
  operations from a Fortran code but using the
  Xwindow protocole

  developed by Jean B. on August 3rd, 1990.

                     VERSION 0.0

*/

/*                global definitions               */

#include <X11/Xlib.h>
#include <X11/Xutil.h>
#include <X11/cursorfont.h>
#include <X11/X.h>
#include <X11/keysym.h>
#include <stdio.h>

#define PIXELS    230
#define MAXDOT    200

Window                 master_window,plot_window,title_window;
Window                 interim_plot_window;
XImage                *bitmap_image;
int                    bitmap_width,bitmap_height;
Window                 input_window,output_window;
Window                 interrupt_window;
Display               *mydisplay;
GC                     master_gc,plot_gc,title_gc;
GC                     interim_plot_gc;
GC                     input_gc,output_gc,icon_gc;
GC                     interrupt_gc;
XEvent                 myevent;
XSizeHints             myhint;
Pixmap                 save_pixmap;
int                    myscreen;
unsigned long          pixels[PIXELS];
Colormap               cmap,screen_cmap;
unsigned int           ncolors;
unsigned long          plane_masks;
Font                   smallfont,mediumfont,largefont,superfont;
XFontStruct           *mediumfont_struct;
int                    mylinestyle;
int                    mylinewidth;
XColor                 color[PIXELS+24],ground_color;
unsigned int           plot_width,plot_height;
unsigned int           depth;                                                
int                    print_count;
FILE                  *fopen(), *fp;
char                   he[512];
char                   output_message[256];
int                    len_output_message;
Cursor                 title_cursor,input_cursor,output_cursor;
Cursor                 interrupt_cursor;
Pixmap                 logo_pixmap;
char                   save_mode;
Cursor                 mycursor,arrowcursor;
unsigned long          myforeground,mybackground,darkbackground,whitebackground;
char                   PostScript;

extern void popen_ ();
extern void pclose_ ();
extern void pcmap_ ();
extern void pplot_ ();
extern void pscale_ ();
extern void psave_ ();
extern void pfill_ ();
extern void pline_ ();
extern void ppen_ ();
extern void psymbol_ ();
extern void pclear_ ();
extern void pthick_ ();
extern void pclip_ ();
extern void pclipoff_();
extern void pcircle_();

static char            dash1[]={3,3};
static char            dash2[]={6,3};
static char            dash3[]={6,6};
static char            dash4[]={9,3};
static char            dash5[]={9,3,3,3};
static int             red[]=  {65535,0,0};
static int             green[]={65535,0,0};
static int             blue[]= {65535,0,0};

static  char  logo[] = {0};

struct scaling { int dummy; 
  int   xpmin; int   xpmax; int   ypmin; int   ypmax;
  float xmin;  float xmax;  float ymin;  float ymax;
  char mapping;
  } scale;

struct currently {
  short xposition; short yposition;
  } current;


/* ------------------------------------------------- */
xopen_(window,display,ndisplay)

/* purpose: to open a window and initialize the plot 
   arguments: window ... a integer*4 array of 4 elements
                         containing the upper left corner
                         x and y coordinates in pixel units
                         and the height and width in pixel 
                         units
              display .. character string containing the 
                         name of a X-display 
                         ("rsesa1:0", for instance)
                         a null character (char(0))) means default
                         terminal attached to running CPU

*/

int                    *window;
char                   *display;
int                    ndisplay; 

{

Bool                   contig;
unsigned int           nplanes;
Status                 result,status;
int                    i;
XEvent                 go_event;
XWMHints               wmhint;
int                    cmap_count;
int                    ndirs;

/* sets PostScript option off */

PostScript='n';

/* initializes display and screen */

mydisplay=XOpenDisplay(display);
myscreen=DefaultScreen(mydisplay);
depth=XDefaultDepth(mydisplay,myscreen);
save_mode='f';
arrowcursor=XCreateFontCursor(mydisplay,XC_tcross);

/* load 3 fonts */

smallfont=XLoadFont(mydisplay,"6x10");
mediumfont=XLoadFont(mydisplay,"6x13");
mediumfont_struct=XLoadQueryFont(mydisplay,"6x13");
largefont=XLoadFont(mydisplay,"9x15");
superfont=XLoadFont(mydisplay,"9x15");

/* defines master window geometry */

myhint.x= window[0]-110;
myhint.y= window[1]-15;
myhint.width= window[2]+115;
myhint.height= window[3]+42;
myhint.flags=PPosition|PSize;

plot_width=window[2];
plot_height=window[3];

/* creates general color map */

if (depth==24) 
{
cmap=XCreateColormap(mydisplay,DefaultRootWindow(mydisplay),
                     DefaultVisual(mydisplay,myscreen),
                     AllocNone);
ground_color.red=65000;
ground_color.green=65335;
ground_color.blue=60000;
ground_color.flags=DoRed | DoGreen | DoBlue;
XAllocColor (mydisplay,cmap, &ground_color);
mybackground=ground_color.pixel;
ground_color.red=0;
ground_color.green=0;
ground_color.blue=0;
ground_color.flags=DoRed | DoGreen | DoBlue;
XAllocColor (mydisplay,cmap, &ground_color);
myforeground=ground_color.pixel;
ground_color.red=50000;
ground_color.green=55000;
ground_color.blue=45000;
ground_color.flags=DoRed | DoGreen | DoBlue;
XAllocColor (mydisplay,cmap, &ground_color);
darkbackground=ground_color.pixel;
ground_color.red=65335;
ground_color.green=65335;
ground_color.blue=65335;
ground_color.flags=DoRed | DoGreen | DoBlue;
XAllocColor (mydisplay,cmap, &ground_color);
whitebackground=ground_color.pixel;
}

else
{
mybackground=230;
myforeground=231;
darkbackground=232;
whitebackground=233;
cmap=XCreateColormap(mydisplay,DefaultRootWindow(mydisplay),
                     DefaultVisual(mydisplay,myscreen),
                     AllocAll);
}

screen_cmap=XDefaultColormap (mydisplay,myscreen);

/* create master window */

master_window=XCreateSimpleWindow(mydisplay,
  DefaultRootWindow (mydisplay),
  myhint.x,myhint.y,myhint.width,myhint.height,
  2,myforeground,darkbackground);

logo_pixmap=XCreatePixmapFromBitmapData (mydisplay,
         master_window,
         logo,100,120,
         myforeground,mybackground,
         depth
         );

wmhint.flags=InputHint;
wmhint.input=True;
XSetWMHints (mydisplay,master_window, &wmhint);

XStoreName (mydisplay,master_window,"xplot\n");
XSetIconName (mydisplay,master_window,"xplot\n");

/* create graphic content */

master_gc=XCreateGC(mydisplay,master_window,0,0);

/* enable input events*/

XSelectInput(mydisplay,master_window,
             ExposureMask|EnterWindowMask|LeaveWindowMask|KeyPressMask);

/* make the window visible */

XMapRaised(mydisplay,master_window);

/* same thing for the bar window */

title_window=XCreateSimpleWindow(mydisplay,
  master_window,
  14,1,myhint.width-18,10,
  1,myforeground,mybackground);

/* create graphic content */

title_gc=XCreateGC(mydisplay,title_window,0,0);

XSetForeground(mydisplay,title_gc,myforeground);

XSelectInput(mydisplay,title_window,
             ButtonPressMask|PointerMotionMask);

/* make the window visible */

XMapRaised(mydisplay,title_window);

/* define special cursor in bar window */

title_cursor=XCreateFontCursor(mydisplay,XC_dot);
XDefineCursor(mydisplay,title_window,title_cursor);

/* same thing for the interrupt window */

interrupt_window=XCreateSimpleWindow(mydisplay,
  master_window,
  2,1,10,10,
  1,myforeground,mybackground);

/* create graphic content */

interrupt_gc=XCreateGC(mydisplay,interrupt_window,0,0);

XSelectInput(mydisplay,interrupt_window,
             ButtonPressMask);

/* make the window visible */

XMapRaised(mydisplay,interrupt_window);

/* define special cursor in bar window */

interrupt_cursor=XCreateFontCursor(mydisplay,XC_pirate);
XDefineCursor(mydisplay,interrupt_window,interrupt_cursor);

/* same thing for the input window */

input_window=XCreateSimpleWindow(mydisplay,
  master_window,
  myhint.width/2+3,myhint.height-22,(myhint.width-15)/2,17,
  1,myforeground,mybackground);

/* create graphic content */

input_gc=XCreateGC(mydisplay,input_window,0,0);

XSelectInput(mydisplay,input_window,
             KeyPressMask);

/* make the window visible */

XMapRaised(mydisplay,input_window);

/* define special cursor in input window */

input_cursor=XCreateFontCursor(mydisplay,XC_left_side);
XDefineCursor(mydisplay,input_window,input_cursor);
XSetForeground(mydisplay,input_gc,myforeground);

/* same thing for the output window */

output_window=XCreateSimpleWindow(mydisplay,
  master_window,
  2,myhint.height-22,(myhint.width-15)/2,17,
  1,myforeground,mybackground);

/* create graphic content */

output_gc=XCreateGC(mydisplay,output_window,0,0);

/* make the window visible */

XMapRaised(mydisplay,output_window);

/* define special cursor in output window */

output_cursor=XCreateFontCursor(mydisplay,XC_gumby);
XDefineCursor(mydisplay,output_window,output_cursor);
XSetForeground(mydisplay,output_gc,myforeground);

/* same thing for the plot window */

/* create window */

plot_window=XCreateSimpleWindow(mydisplay,
  master_window,
  110,15,plot_width,plot_height,
  1,myforeground,whitebackground);

/* create graphic content */

plot_gc=XCreateGC(mydisplay,plot_window,0,0);

/* enable input events*/

XSelectInput(mydisplay,plot_window,
  ButtonPressMask|ExposureMask|
  EnterWindowMask|LeaveWindowMask|PointerMotionMask);

/* make the window visible */

XMapRaised(mydisplay,plot_window);

/* attach color map to plot window */

ncolors=PIXELS;

XSetWindowColormap(mydisplay,master_window,cmap);
for (i=0;i<ncolors;i++) pixels[i]=20+i;

/* loads the fonts */

XSetFont(mydisplay,output_gc,mediumfont);
XSetFont(mydisplay,input_gc,mediumfont);

XSetFillRule(mydisplay,plot_gc,WindingRule);

/* create a pixmap (the size of the window) 
   to save the window content when xsave is called */

save_pixmap=XCreatePixmap
  (mydisplay,plot_window,plot_width,plot_height,depth);

mylinestyle=LineSolid;

scale.mapping='w';

current.xposition=0;
current.yposition=0;

len_output_message=0;

bitmap_width=0;

while (True)
  {
  XNextEvent (mydisplay,&go_event);
  if ( go_event.type = Expose ) return;
  }

}

/*-------------------------------------------------- */
xpson_(fnme,nlen,nfnme)

/* to set PostScript option on 
     when the postcript option is on, all graphic operations
     are also sent to a Postscript file

     fnme: name of the file
     nlen: length of fnme
     nfnme: dummy argument (not required in Fortran)
*/

char     *fnme;
int      *nfnme;
int      *nlen;

{

int    window[2],size_of_window;

window[0]=plot_width; window[1]=plot_height;
size_of_window=2;
popen_ (window,size_of_window,fnme,nlen,nfnme);
PostScript='y';
XSetForeground(mydisplay,interrupt_gc,mybackground);
XSetBackground(mydisplay,interrupt_gc,myforeground);
XFillRectangle (mydisplay,interrupt_window,interrupt_gc,0,0,10,10);

}

/*-------------------------------------------------- */
xpsoff_()

/* to set PostScript option off */

{

pclose_ ();
PostScript='n';
XSetForeground(mydisplay,interrupt_gc,myforeground);
XSetBackground(mydisplay,interrupt_gc,mybackground);
XFillRectangle (mydisplay,interrupt_window,interrupt_gc,0,0,10,10);

}

/* ------------------------------------------------- */
xpspause_()

/* to pause PostScript output to PostScript file*/

{

PostScript='n';

}

/* ------------------------------------------------- */
xpsresume_()

/* to resume PostScript output to PstScript file */

{

PostScript='y';

}

/* ------------------------------------------------- */
xclose_()
  
/* to close the plot and erase the window */

{          
/* release the save pixmap */

XFreePixmap(mydisplay,save_pixmap);
XFreePixmap(mydisplay,logo_pixmap);

/* free the window and destroy its mapping on screen */

/* XFreeColormap(mydisplay,cmap); */
XFreeGC(mydisplay,master_gc);
XDestroyWindow(mydisplay,master_window);
XCloseDisplay(mydisplay);

}
       
/* ------------------------------------------------- */
xmode_(mode,lmode)

/* purpose: to select the save mode
   argument: mode ... f = fast mode
                      s = save mode
*/

char    *mode;
int     *lmode;

{

save_mode= *mode;

}

/* ------------------------------------------------- */
xpause_(xcursor,ycursor,button)

/* purpose: to pause plotting and follow cursor position
   arguments: xcursor ... x-position of cursor when
                          mouse button was pressed
              ycursor ... y-position of cursor when
                          mouse button was pressed
              button .... mouse button number
                          (neg.=button number)
                          (pos.=ascii code of key pressed)
*/

float       *xcursor, *ycursor;
int         *button;

{                      

int     wait;
int     xm,ym;
GC      mask_gc; 

xm=0;ym=0;

mask_gc=XCreateGC (mydisplay,plot_window,0,0);
XSetFunction (mydisplay,mask_gc,GXinvert);
XSetPlaneMask (mydisplay,mask_gc,AllPlanes);

XDrawLine(mydisplay,plot_window,mask_gc,
          xm,0,xm,plot_height);
XDrawLine(mydisplay,plot_window,mask_gc,
          0,ym,plot_width,ym);

XDefineCursor(mydisplay,plot_window,arrowcursor);

wait=0; 
while (wait==0)
  {
/* wait for an event */

  XNextEvent (mydisplay,&myevent);

  switch (myevent.type){

/* the window is brought back to the surface */
  case Expose:
    if (myevent.xexpose.window == plot_window)
    XCopyArea(mydisplay,save_pixmap,
              plot_window,plot_gc,
              myevent.xexpose.x,myevent.xexpose.y,
              myevent.xexpose.width,myevent.xexpose.height,
              myevent.xexpose.x,myevent.xexpose.y);
    if (myevent.xexpose.window == output_window)
    {
    XClearWindow(mydisplay,output_window);
    XDrawString(mydisplay,output_window,output_gc,
                2,20,output_message,len_output_message);
    }
    break;

  case MappingNotify:
    XRefreshKeyboardMapping( &myevent);
    break;

/* during cursor movement */
  case MotionNotify:
    XDrawLine(mydisplay,plot_window,mask_gc,
              xm,0,xm,plot_height);
    XDrawLine(mydisplay,plot_window,mask_gc,
              0,ym,plot_width,ym);
    xm=myevent.xmotion.x; ym=myevent.xmotion.y;
    XDrawLine(mydisplay,plot_window,mask_gc,
              xm,0,xm,plot_height);
    XDrawLine(mydisplay,plot_window,mask_gc,
              0,ym,plot_width,ym);
    break;

/* a mouse button is pressed */
  case ButtonPress:
    switch (scale.mapping) {
    case ('w'):
      *xcursor = myevent.xbutton.x; *ycursor = myevent.xbutton.y;
      break;
    case ('u'):
      *xcursor=((myevent.xbutton.x-scale.xpmin)
               *(scale.xmax-scale.xmin))
               /(scale.xpmax-scale.xpmin)+scale.xmin;
      *ycursor=((myevent.xbutton.y-scale.ypmin)
               *(scale.ymax-scale.ymin))
               /(scale.ypmax-scale.ypmin)+scale.ymin;
      break;
      }                     
    *button= myevent.xbutton.button;
    XDrawLine(mydisplay,plot_window,mask_gc,
              xm,0,xm,plot_height);
    XDrawLine(mydisplay,plot_window,mask_gc,
              0,ym,plot_width,ym);
    wait=1;
    break;

    }

  }

XDefineCursor(mydisplay,plot_window,mycursor);

}

/* ------------------------------------------------- */
xscale_(xpmi,xpma,ypmi,ypma,xmi,xma,ymi,yma,map)

/* purpose: to define a user system of coordinate and 
            swith between pixel coordinates and user
            coordinates 
   arguments: xpmi,...,ypma ..... mapped area in pixel units
              xmi,...yma ........ mapped are in user units
              map ............... =1 user mapping
                                  =0 pixel mapping
*/

int       *xpmi,*xpma,*ypmi,*ypma;
float     *xmi,*xma,*ymi,*yma;
int       *map;

{

if (PostScript == 'y') pscale_ (xpmi,xpma,ypmi,ypma,xmi,xma,ymi,yma,map);

  switch (*map) {
 
  case (0):
    scale.mapping='w';
    break;
  case (1):
    scale.mapping='u';
    scale.xmin = *xmi;
    scale.xmax = *xma;
    scale.ymin = *ymi;
    scale.ymax = *yma;
    scale.xpmin = *xpmi;
    scale.xpmax = *xpma;
    scale.ypmin = *ypmi;
    scale.ypmax = *ypma;
    break;

  }

}

/* ------------------------------------------------- */
xplot_(x,y,pen)

/* purpose: to draw a line from current position to
            new position (x,y)
   arguments: x ...... new x-position
              y ...... new y-position
              pen .... 2=solid line
                       3=pen up
                       4 to 8=use dash patterns
*/

float     *x,*y;
int       *pen;

{

short           xplot, yplot;
float           ratio,dx,dy;

if (PostScript == 'y') pplot_(x,y,pen);

/* calculate xplot,yplot depending on the mapping
   (user or pixel) */

  switch (scale.mapping) {

  case ('w'):
    xplot = *x; yplot = *y;
    break;
 
  case ('u'):
    xplot=( *x-scale.xmin)/(scale.xmax-scale.xmin)
         *(scale.xpmax-scale.xpmin)+scale.xpmin;
    yplot=( *y-scale.ymin)/(scale.ymax-scale.ymin)
         *(scale.ypmax-scale.ypmin)+scale.ypmin;
    break;

  }
  switch (*pen) {

/* dash pattern 5 */

  case (8):
    XSetDashes (mydisplay,plot_gc,0,dash5,4);
    mylinestyle=LineOnOffDash;
    XSetLineAttributes (mydisplay,plot_gc,mylinewidth,mylinestyle,
                        CapProjecting,JoinBevel);
    break;

/* dash pattern 4 */

  case (7):
    XSetDashes (mydisplay,plot_gc,0,dash4,2);
    mylinestyle=LineOnOffDash;
    XSetLineAttributes (mydisplay,plot_gc,mylinewidth,mylinestyle,
                        CapProjecting,JoinBevel);
    break;

/* dash pattern 3 */

  case (6):
    XSetDashes (mydisplay,plot_gc,0,dash3,2);
    mylinestyle=LineOnOffDash;
    XSetLineAttributes (mydisplay,plot_gc,mylinewidth,mylinestyle,
                        CapProjecting,JoinBevel);
    break;

/* dash pattern 2 */

  case (5):
    XSetDashes (mydisplay,plot_gc,0,dash2,2);
    mylinestyle=LineOnOffDash;
    XSetLineAttributes (mydisplay,plot_gc,mylinewidth,mylinestyle,
                        CapProjecting,JoinBevel);
    break;

/* dash pattern 1 */

  case (4):
    XSetDashes (mydisplay,plot_gc,0,dash1,2);
    mylinestyle=LineOnOffDash;
    XSetLineAttributes (mydisplay,plot_gc,mylinewidth,mylinestyle,
                        CapProjecting,JoinBevel);
    break;

/* pen up */

  case (3):
    break;

/* solid line */

  case (2):
    mylinestyle=LineSolid;
    XSetLineAttributes (mydisplay,plot_gc,mylinewidth,mylinestyle,
                        CapProjecting,JoinBevel);
    break;

  }

  if ( *pen != 3)
  XDrawLine (mydisplay,plot_window,plot_gc,
             current.xposition,current.yposition,
             xplot,yplot);

  if ( *pen != 3 && save_mode == 's')
  XDrawLine (mydisplay,save_pixmap,plot_gc,
             current.xposition,current.yposition,
             xplot,yplot);

/* update the current position */

  current.xposition=xplot; current.yposition=yplot;

}

/* ------------------------------------------------- */
xsave_()

/* purpose: to save the present window graphic content
            to be brought back to the surface when the
            window is unobscured 
*/

{

if (PostScript == 'y') psave_ ();

if (bitmap_width != 0) return;  

XFlush(mydisplay);
                                                                
if (save_mode == 'f') 
  XCopyArea
  (mydisplay,plot_window,save_pixmap,plot_gc,0,0,
   plot_width,plot_height,0,0);

}

/* ------------------------------------------------- */
xfill_(x,y,n)

/* purpose: to fill a (x,y) polygon 
   arguments: x,y ...... coordinate arrays
              n ........ number of points in x and y

        *note: maximum of 100000 points allowed*

*/

float     *x, *y;
int       *n;

{

XPoint    figure[100000];
int       i,npoint;

if (PostScript == 'y') pfill_ (x,y,n);

npoint = *n;                  
              
if ( npoint > 100000) return;

/* calculate figure depending on the mapping
   (user or pixel) */

switch (scale.mapping) {

  case ('w'):
    for (i = 0 ; i < npoint ; i++) 
      {
      figure[i].x = x[i]; figure[i].y = y[i];
      }
    break;
 
  case ('u'):
    for (i = 0 ; i < npoint ; i++) 
      {      
      figure[i].x=( x[i]-scale.xmin)/(scale.xmax-scale.xmin)
           *(scale.xpmax-scale.xpmin)+scale.xpmin;
      figure[i].y=( y[i]-scale.ymin)/(scale.ymax-scale.ymin)
           *(scale.ypmax-scale.ypmin)+scale.ypmin;
      }
    break;

  }

XFillPolygon (mydisplay,plot_window,plot_gc,figure,npoint,
              Complex,CoordModeOrigin);

if (save_mode == 's') 
XFillPolygon (mydisplay,save_pixmap,plot_gc,figure,npoint,
              Complex,CoordModeOrigin);

}

/* ------------------------------------------------- */
xline_(x,y,n,pen)

/* purpose: to draw a line passing trough a series
            of points (this subroutine should be used
            rather than xplot to connect more than
            two points with a line, not only because 
            it is faster but also because it respects
            dash patterns across segments
   arguments: x,y ...... corrdinates of points
              n ........ number of points
              pen ...... identical to pen in xplot

        *note: maximum of 100000 points allowed*

*/

float     *x, *y;
int       *n, *pen;

{

XPoint    figure[100000];
int       i,npoint;

if (PostScript == 'y') pline_ (x,y,n,pen);

npoint = *n;                  
              
if ( npoint > 100000) return;

/* calculate figure depending on the mapping
   (user or pixel) and update current pen position*/

switch (scale.mapping) {

  case ('w'):
    for (i = 0 ; i < npoint ; i++) 
      {
      figure[i].x = x[i]; figure[i].y = y[i];
      }
    current.xposition=figure[npoint-1].x;
    current.yposition=figure[npoint-1].y;
    break;
 
  case ('u'):
    for (i = 0 ; i < npoint ; i++) 
      {      
      figure[i].x=( x[i]-scale.xmin)/(scale.xmax-scale.xmin)
           *(scale.xpmax-scale.xpmin)+scale.xpmin;
      figure[i].y=( y[i]-scale.ymin)/(scale.ymax-scale.ymin)
           *(scale.ypmax-scale.ypmin)+scale.ypmin;
      }
    current.xposition=figure[npoint-1].x;
    current.yposition=figure[npoint-1].y;
    break;

  }

  switch (*pen) {

/* dash pattern 5 */

  case (8):
    XSetDashes (mydisplay,plot_gc,0,dash5,4);
    mylinestyle=LineOnOffDash;
    XSetLineAttributes (mydisplay,plot_gc,mylinewidth,mylinestyle,
                        CapProjecting,JoinBevel);
    break;

/* dash pattern 4 */

  case (7):
    XSetDashes (mydisplay,plot_gc,0,dash4,2);
    mylinestyle=LineOnOffDash;
    XSetLineAttributes (mydisplay,plot_gc,mylinewidth,mylinestyle,
                        CapProjecting,JoinBevel);
    break;

/* dash pattern 3 */

  case (6):
    XSetDashes (mydisplay,plot_gc,0,dash3,2);
    mylinestyle=LineOnOffDash;
    XSetLineAttributes (mydisplay,plot_gc,mylinewidth,mylinestyle,
                        CapProjecting,JoinBevel);
    break;

/* dash pattern 2 */

  case (5):
    XSetDashes (mydisplay,plot_gc,0,dash2,2);
    mylinestyle=LineOnOffDash;
    XSetLineAttributes (mydisplay,plot_gc,mylinewidth,mylinestyle,
                        CapProjecting,JoinBevel);
    break;

/* dash pattern 1 */

  case (4):
    XSetDashes (mydisplay,plot_gc,0,dash1,2);
    mylinestyle=LineOnOffDash;
    XSetLineAttributes (mydisplay,plot_gc,mylinewidth,mylinestyle,
                        CapProjecting,JoinBevel);
    break;

/* pen up */

  case (3):
    break;

/* solid line */

  case (2):
    mylinestyle=LineSolid;
    XSetLineAttributes (mydisplay,plot_gc,mylinewidth,mylinestyle,
                        CapProjecting,JoinBevel);
    break;

  }

if ( *pen != 3)
XDrawLines (mydisplay,plot_window,plot_gc,figure,npoint,
            CoordModeOrigin);

if ( *pen != 3 && save_mode == 's')
XDrawLines (mydisplay,save_pixmap,plot_gc,figure,npoint,
            CoordModeOrigin);

}

/* ------------------------------------------------- */
xcmap_(r,g,b)

/* purpose: to define a new color map; one must define 230 color
            at a time (r,g,b must be dimensioned (1:230) in
            FORTRAN and all values in r,g and b must be defined
            minimum intensity:0 ; maximum intensity: 255.
   arguments: r ...... 230 red values
              g ...... 230 green values
              b ...... 230 blue values
*/

int      *r, *g, *b;               
                    
{

int                i;
int                ncolmin,ncolmax;
XColor             color_get;

if (PostScript == 'y') pcmap_ (r,g,b);
                  
ncolmin=0;                                            
ncolmax=ncolors;

   for (i=0;i<20;i++)
  {
  color_get.pixel=i;
  XQueryColor (mydisplay,screen_cmap, &color_get);
  color[ncolmax+i+1].red=color_get.red;
  color[ncolmax+i+1].green=color_get.green;
  color[ncolmax+i+1].blue=color_get.blue;
  color[ncolmax+i+1].flags= DoRed | DoGreen | DoBlue;
  color[ncolmax+i+1].pixel=i;
  }

  for (i=ncolmin;i<ncolmax;i++)
  {
  color[i].red  =  (65535*r[i])/255;
  color[i].green=  (65535*g[i])/255;
  color[i].blue =  (65535*b[i])/255;
  color[i].flags=  DoRed | DoGreen | DoBlue;
  color[i].pixel=  pixels[i];
  }

color[ncolmax+20].red=65000;
color[ncolmax+20].green=65335;
color[ncolmax+20].blue=60000;
color[ncolmax+20].flags=DoRed | DoGreen | DoBlue;
color[ncolmax+20].pixel=ncolmax;

color[ncolmax+21].red=0;
color[ncolmax+21].green=0;
color[ncolmax+21].blue=0;
color[ncolmax+21].flags=DoRed | DoGreen | DoBlue;
color[ncolmax+21].pixel=ncolmax+1;

color[ncolmax+22].red=50000;
color[ncolmax+22].green=55000;
color[ncolmax+22].blue=45000;
color[ncolmax+22].flags=DoRed | DoGreen | DoBlue;
color[ncolmax+22].pixel=ncolmax+2;

color[ncolmax+23].red=65535;
color[ncolmax+23].green=65535;
color[ncolmax+23].blue=65535;
color[ncolmax+23].flags=DoRed | DoGreen | DoBlue;
color[ncolmax+23].pixel=ncolmax+3;

if (depth==24)
{
for (i=ncolmin;i<ncolmax;i++) {color[i].pixel=0;
                                 XAllocColor (mydisplay,cmap, &color[i]);}
for (i=ncolmin;i<ncolmax;i++) pixels[i]=color[i].pixel;
}
else
XStoreColors (mydisplay,cmap,color,ncolors+24);

XFlush(mydisplay);

}

/* ------------------------------------------------- */
xpen_(pen)

/* purpose: to change current pen (color)
   argument: pen ..... new pen (color)
*/

int    *pen;

{

int    color;

if (PostScript == 'y') ppen_ (pen);

color = *pen - 1;

if (color > ncolors || color < 0) return;

XSetForeground (mydisplay,plot_gc,pixels[color]);

}

/* ------------------------------------------------- */
xcursor_()

/* purpose: to activate the cursor */

{

mycursor=XCreateFontCursor(mydisplay,XC_cross_reverse);

XDefineCursor(mydisplay,plot_window,mycursor);

}

/* ------------------------------------------------- */
xsymbol_(x,y,size,string,nlen,dummy)

/* purpose: to draw a string 
   arguments: x,y ..... coordinate of lower left corner of
                        string
              size .... size of string (1, 2 or 3)
                        it will select pre-defined fonts
              nlen .... length of string
            
*/

float      *x, *y;
char       *string;
int        *nlen;
int         dummy;
int        *size;

{
                                    
short    xstring,ystring;

if (PostScript == 'y') psymbol_ (x,y,size,string,nlen,dummy);
 
/* pick up the font according to size */

switch ( *size) 
  {
  case (1):
    XSetFont(mydisplay,plot_gc,smallfont);
    break;

  case (2):
    XSetFont(mydisplay,plot_gc,mediumfont);
    break;

  case (3):
    XSetFont(mydisplay,plot_gc,largefont);
    break;

  case (4):
    XSetFont(mydisplay,plot_gc,superfont);
    break;                       
  }

/* find position of string depending on mapping */

switch (scale.mapping) {

  case ('w'):
    xstring = *x; ystring = *y;
    break;
 
  case ('u'):
    xstring=( *x-scale.xmin)/(scale.xmax-scale.xmin)
         *(scale.xpmax-scale.xpmin)+scale.xpmin;
    ystring=( *y-scale.ymin)/(scale.ymax-scale.ymin)
         *(scale.ypmax-scale.ypmin)+scale.ypmin;
    break;

  }

/* draw the string */

XDrawString(mydisplay,plot_window,plot_gc,
            xstring,ystring,string,*nlen);

if (save_mode == 's') 
XDrawString(mydisplay,save_pixmap,plot_gc,
            xstring,ystring,string,*nlen);

}

/* ------------------------------------------------- */
xclear_()

/* purpose: to clear the window */

{

if (PostScript == 'y') pclear_ ();

XClearWindow(mydisplay,plot_window);
                                         
  if (save_mode == 's') 
  {
  XSetForeground(mydisplay,plot_gc,mybackground);
  XFillRectangle(mydisplay,save_pixmap,
       plot_gc,0,0,plot_width,plot_height);
  }

}

/* ------------------------------------------------- */
xthick_(width)

/* purpose: to change line width
   argument: width ... thickness in pixel

*/

int      *width;

{

if (PostScript == 'y') pthick_ (width);

mylinewidth= *width;
XSetLineAttributes (mydisplay,plot_gc,mylinewidth,
                    mylinestyle,CapProjecting,JoinBevel);

}

/* ------------------------------------------------- */
xwhere_(x,y)

/* purpose: returns current pen position 
   arguments: x,y .... current pen coordinates
*/

float  *x, *y;

{

switch (scale.mapping)

  {

  case ('w'):
    *x=current.xposition; *y=current.yposition;
    break;
 
  case ('u'):
    *x=(current.xposition-scale.xpmin)/(scale.xpmax-scale.xpmin)*
       (scale.xmax-scale.xmin)+scale.xmin;
    *y=(current.yposition-scale.ypmin)/(scale.ypmax-scale.ypmin)*
       (scale.ymax-scale.ymin)+scale.ymin;

    break;
   
  }

}

/* ------------------------------------------------- */
xpopup_()                                     

/* purpose: to bring the window back to the surface */

{

if (bitmap_width != 0) return;

XMapRaised (mydisplay,master_window);
XCopyArea(mydisplay,save_pixmap,
          plot_window,plot_gc,
          0,0,plot_width,plot_height,0,0);

}

/* ------------------------------------------------- */
xpoints_(x,y,n)

/* purpose to paint a series of pixel

  arguments: x and y are two integer*2 arrays containing the
             coordinates of the pixels to be painted
             n number of points
  note that this function can only be used if the x-y coordinates
  are pixel valued (this routine does not use user coordinates)
*/

short     *x, *y;
int       *n;

{
                                          
XPoint    points[10000];
int       i;

  for ( i=0 ; i< *n ; i++) 
  {
  points[i].x= x[i]; points[i].y= y[i];
  }

XDrawPoints (mydisplay,plot_window,plot_gc,
             points, *n, CoordModeOrigin);

if (save_mode == 's')
XDrawPoints (mydisplay,save_pixmap,plot_gc,
             points, *n, CoordModeOrigin);

}

/* -------------------------------------------------- */
xclip_(xmin_clip,xmax_clip,ymin_clip,ymax_clip)

/* purpose: to define a clipping window
   arguments: xmin_clip ..... minimum x value
              xmax_clip ..... maximum x value
              ymin_clip ..... minimum y value
              ymax_clip ..... maximum y value
         any graphic operation outside this window
         is not performed

*/

float      *xmin_clip, *xmax_clip, *ymin_clip, *ymax_clip;

{

XRectangle rect[1];

if (PostScript == 'y') pclip_ (xmin_clip,xmax_clip,ymin_clip,ymax_clip);

  switch (scale.mapping) {

  case ('w'):
    rect[0].x = *xmin_clip ; rect[0].y = *ymin_clip;
    rect[0].height = *xmax_clip - *xmin_clip;
    rect[0].width = *ymax_clip - *ymin_clip;
    break;
 
  case ('u'):
    rect[0].x=( *xmin_clip-scale.xmin)/(scale.xmax-scale.xmin)
          *(scale.xpmax-scale.xpmin)+scale.xpmin;
    rect[0].y=( *ymin_clip-scale.ymin)/(scale.ymax-scale.ymin)
         *(scale.ypmax-scale.ypmin)+scale.ypmin;
    rect[0].width=( *xmax_clip- *xmin_clip)/(scale.xmax-scale.xmin)
          *(scale.xpmax-scale.xpmin);
    rect[0].height=( *ymax_clip- *ymin_clip)/(scale.ymax-scale.ymin)
         *(scale.ypmax-scale.ypmin);
    break;

  }

XSetClipRectangles (mydisplay,plot_gc,0,0,
                    rect,1,Unsorted);

}

/* -------------------------------------------------- */
xclipoff_()

/* purpose: to reset the clipping window to its default:
            the entire window */

{

if (PostScript == 'y') pclipoff_ ();

XSetClipMask (mydisplay,plot_gc,None);

}
/* -------------------------------------------------- */
int xchoice_(choice,nchoice,dummy)

/* purpose: to interactively make a choice 
            xchoice is a integer function which return
            the choice made by the user between a series
            of options displayed in small windows on the
            left side of the main window 
   arguments: choice ..... array of 20 long character string
                           describing the various choices
              nchoice .... number of possible choices 
                           (dimension of choice)

      * nchoice < 20  *

*/

int     *nchoice;
char    *choice;
int     dummy;
 
#define  CHOICELENGTH  20
#define  MAXCHOICE     20
                 
{
 
Window        subwindow[MAXCHOICE];
GC            subgc[MAXCHOICE];
XFontStruct  *mediumfontstruct;
int           i,j;
int           width,maxwidth;
int           height;
XEvent        subevent;
XEvent        created;
char          string[CHOICELENGTH];
short         xstring,ystring;
Cursor        subcursor;
int           nl;
Window        focus;
int           revert_to;
Bool          check;

/* create a new cursor and define a font */

subcursor=XCreateFontCursor(mydisplay,XC_hand2);
nl=CHOICELENGTH;

height=30;
if (plot_height/ *nchoice < height) height=plot_height/ *nchoice;

/* create the small windows */

  for ( i = 0 ; i < *nchoice ; i++ )
  {
  subwindow[i]=XCreateSimpleWindow(mydisplay,
               master_window,
               3,15+height*i+2,100,height-4,
               1,myforeground,mybackground);

  subgc[i]=XCreateGC(mydisplay,subwindow[i],0,0);

  XSelectInput(mydisplay,subwindow[i],
               ButtonPressMask|ExposureMask|
               EnterWindowMask|LeaveWindowMask);

  XMapRaised(mydisplay,subwindow[i]);

  XSetFont(mydisplay,subgc[i],smallfont);
  XDefineCursor(mydisplay,subwindow[i],subcursor);
  XSetForeground(mydisplay,subgc[i],myforeground);

  check=True;
  while (check)
    {
    XNextEvent (mydisplay, &created);
    if ( created.type = Expose || created.xexpose.window == subwindow[i])
     check=False;
    }

    for ( j=0 ; j<nl ; j++ ) string[j]=choice[ i * nl + j ];
    string[nl]='\0';
    xstring=2; ystring=height*2/3;
    XDrawString(mydisplay,subwindow[i],subgc[i],
                xstring,ystring,string,nl);

  }

/* wait for an event */

while (True)
  {
  XNextEvent (mydisplay,&subevent);

  switch (subevent.type)
    {
/* a choice is made */

    case ButtonPress:
      for (i=0;i< *nchoice;i++)
        if (subevent.xbutton.window == subwindow[i])        
          {                                    
          for (j=0 ; j< *nchoice; j++)
            {
            XFreeGC(mydisplay,subgc[j]);
            XDestroyWindow(mydisplay,subwindow[j]);
            }
          XFlush;
          return(i+1);
          }
      break;

    case MappingNotify:
      XRefreshKeyboardMapping( &subevent);
      break;

/* switch black and white if the cursor enters the window */

    case EnterNotify:                 
      for (i=0;i< *nchoice;i++)
        {
        if (subevent.xcrossing.window == subwindow[i])        
          {
          XFillRectangle(mydisplay,subwindow[i],subgc[i],
                         0,0,100,height);
          XSetForeground(mydisplay,subgc[i],mybackground);
          for ( j=0 ; j<nl ; j++ ) string[j]=choice[ i * nl + j ];
          XDrawString(mydisplay,subwindow[i],subgc[i],
                      xstring,ystring,string,nl);
          }
        }
      break;

/* switch black and white if the cursor leaves the window */

    case LeaveNotify:
      for (i=0;i< *nchoice;i++)
        {
        if (subevent.xcrossing.window == subwindow[i])        
          {
          XFillRectangle(mydisplay,subwindow[i],subgc[i],
                         0,0,100,height);
          XSetForeground(mydisplay,subgc[i],myforeground);
          for ( j=0 ; j<nl ; j++ ) string[j]=choice[ i * nl + j ];
          XDrawString(mydisplay,subwindow[i],subgc[i],
                      xstring,ystring,string,nl);
          }
        }
      break;

/* redraws in case of reexposure */

    case Expose:
      if (subevent.xexpose.window == output_window)
      {
      XClearWindow(mydisplay,output_window);
      XDrawString(mydisplay,output_window,output_gc,
                  2,20,output_message,len_output_message);
      }
      if (subevent.xexpose.window == plot_window)
      XCopyArea(mydisplay,save_pixmap,
                plot_window,plot_gc,
                subevent.xexpose.x,subevent.xexpose.y,
                subevent.xexpose.width,subevent.xexpose.height,
                subevent.xexpose.x,subevent.xexpose.y);
      for (i=0;i< *nchoice;i++)
        {
        for ( j=0 ; j<nl ; j++ ) string[j]=choice[ i * nl + j ];
        XSetForeground(mydisplay,subgc[i],myforeground);
        XDrawString(mydisplay,subwindow[i],subgc[i],
                      xstring,ystring,string,nl);
        }
      break;
   }
  }

}

/* -------------------------------------------------- */
xbar_(ratio)
                         
/* subroutine to move the cursor in the bar window
   and return its position to the program */

float         *ratio;

{                            

XEvent        event;
int           status;

while (True)
  {
  XNextEvent (mydisplay, &event);

  switch (event.type)

    {

    case MotionNotify:

      if (event.xmotion.window == title_window)
        {                                        
        XClearWindow(mydisplay,title_window);
        XFillRectangle(mydisplay,title_window,title_gc,
                       event.xmotion.x-5,0,10,10);
        *ratio=event.xmotion.x;
        *ratio= *ratio/(myhint.width-18);
        return;     
        }

      break;

    case ButtonPress:
                   
      *ratio= -1.;
      return;

    }

  }

}

/* -------------------------------------------------- */
xtrack_(xpos,ypos,flag)

/* subroutine to track the position of the cursor
   in the plot window 
   arguments: xpos and ypos are the position of the cursor
              flag is always -1 unless a button is pushed
                   in which case the button number is returned in
                   flag
*/
                         
float         *xpos, *ypos;
int           *flag;

{                            

XEvent        event;

while (True)
  {
  XNextEvent (mydisplay, &event);

  switch (event.type)

    {

      case Expose:
      if (event.xexpose.window == plot_window)
      XCopyArea(mydisplay,save_pixmap,
                plot_window,plot_gc,
                event.xexpose.x,event.xexpose.y,
                event.xexpose.width,event.xexpose.height,
                event.xexpose.x,event.xexpose.y);
      if (event.xexpose.window == output_window)
      {
      XClearWindow(mydisplay,output_window);
      XDrawString(mydisplay,output_window,output_gc,
                  2,20,output_message,len_output_message);
      }
      break;

    case MappingNotify:
      XRefreshKeyboardMapping( &event);
      break;

    case MotionNotify:
      if (event.xmotion.window == plot_window)
        {                                        
        switch (scale.mapping) {
        case ('w'):
          *xpos = event.xmotion.x; *ypos = event.xmotion.y;
          break;
        case ('u'):
          *xpos=((event.xmotion.x-scale.xpmin)
                *(scale.xmax-scale.xmin))
                /(scale.xpmax-scale.xpmin)+scale.xmin;
          *ypos=((event.xmotion.y-scale.ypmin)
                *(scale.ymax-scale.ymin))
                 /(scale.ypmax-scale.ypmin)+scale.ymin;
          break;
          }                     
        *flag= -1;
        return;     
        }
      break;

    case ButtonPress:
        *flag=event.xbutton.button;
      return;

    }

  }

}

/* -------------------------------------------------- */
xprint_(fnme,nlen,nfnme)

/* subroutine to dump the contain of the window onto a pixel map in
   Postscript format (to be printed later)

   arguments: fnme: file name
              nfnme: length of fnme
              nfnme: dummy argument
*/

char     *fnme;
int      *nfnme;
int      *nlen;

{
                               
XImage    *image;
int        i,j,nl;
char       file[256];
XColor     print_color[800];
int        grey[800];
char       he2[2];

xpopup_();   

image=XGetImage (mydisplay,save_pixmap,0,0,
                 plot_width,plot_height,AllPlanes,ZPixmap);

for ( i=0 ; i< *nlen ; i++) file[i]=fnme[i];
file[ *nlen]='\0';
fp=fopen(file,"w");
  if (fp==NULL) 
  {
  printf("Could not open file %s\n",file);
  return;
  }

print_count=0;

  for ( i=0 ; i<256 ; i++) 
  {
  sprintf(he2,"%2x",i);
  he[i*2]=he2[0];he[i*2+1]=he2[1];
  }
for ( i=0 ; i<16 ; i++) he[i*2]='0';

fprintf(fp,"%%!2DMap\n");
fprintf(fp,"%%%%%e%e%e%e\n",
      (378-(plot_width*378)/1024)*2.54/72.,
      (378-(plot_width*378)/1024+(plot_width*756)/1024)*2.54/72.,
      (288-(plot_height*288)/800)*2.54/72.,
      (288-(plot_height*288)/800+(plot_height*576)/800)*2.54/72.);
fprintf(fp,"gsave\n");
fprintf(fp,"/picstr %d string def\n",plot_height);
fprintf(fp,"%d %d translate\n",288-(plot_height*288)/800,
                               378-(plot_width*378)/1024);
fprintf(fp,"%d %d scale\n",(plot_height*576)/800,
                           (plot_width*756/1024));
fprintf(fp,"%d %d %d\n",plot_height,plot_width,8);
fprintf(fp,"[%d 0 0 %d 0 %d]\n",plot_height,-plot_width,
                                plot_width);
fprintf(fp,"{currentfile picstr readhexstring pop}\n");
fprintf(fp,"image\n");

  for ( i = plot_width-1 ; i > -1 ; i-- )
  {
    if (xinterrupt_()==1) 
    {
    fclose(fp);
    return;
    }
  XFillRectangle(mydisplay,title_window,title_gc,
                 98+i-1,0,2,20);
    for (j=0;j<plot_height;j++) 
    print_color[j].pixel=XGetPixel(image,i,j);
  XQueryColors (mydisplay,cmap,print_color,plot_height);
    for (j=0;j<plot_height;j++)
  grey[j]=
   ( .35*print_color[j].red+
     .55*print_color[j].green+
     .10*print_color[j].blue)/256;
  pp(grey,plot_height);
  }

fprintf(fp,"\nstroke\n");
fprintf(fp,"showpage\n");
fprintf(fp,"grestore\n");
fclose(fp);

XDestroyImage (image);
XClearWindow(mydisplay,title_window);

}

/* -------------------------------------------------- */
xcolorprint_(fnme,nlen,nfnme)

/* subroutine to dump the contain of the window onto a pixel map in
   Postscript format (to be printed later)

   arguments: fnme: file name
              nfnme: length of fnme
              nfnme: dummy argument
   same as xprint but produces a color image
*/
char     *fnme;
int      *nfnme;
int      *nlen;

{
                               
XImage    *image;
int        i,j,nl;
char       file[256];
XColor     print_color[800];
int        grey[800];
char       he2[2];

xpopup_();   

image=XGetImage (mydisplay,save_pixmap,0,0,
                 plot_width,plot_height,AllPlanes,ZPixmap);

for ( i=0 ; i< *nlen ; i++) file[i]=fnme[i];
file[ *nlen]='\0';
fp=fopen(file,"w");
  if (fp==NULL) 
  {
  printf("Could not open file %s\n",file);
  return;
  }

print_count=0;

  for ( i=0 ; i<256 ; i++) 
  {
  sprintf(he2,"%2x",i);
  he[i*2]=he2[0];he[i*2+1]=he2[1];
  }
for ( i=0 ; i<16 ; i++) he[i*2]='0';

fprintf(fp,"%%!2DMap\n");
fprintf(fp,"%%%%%e%e%e%e\n",
      (378-(plot_width*378)/1024)*2.54/72.,
      (378-(plot_width*378)/1024+(plot_width*756)/1024)*2.54/72.,
      (288-(plot_height*288)/800)*2.54/72.,
      (288-(plot_height*288)/800+(plot_height*576)/800)*2.54/72.);
fprintf(fp,"gsave\n");
fprintf(fp,"/redstr %d string def\n",plot_height);
fprintf(fp,"/grestr %d string def\n",plot_height);
fprintf(fp,"/blustr %d string def\n",plot_height);
fprintf(fp,"%d %d translate\n",288-(plot_height*288)/800,
                               378-(plot_width*378)/1024);
fprintf(fp,"%d %d scale\n",(plot_height*576)/800,
                           (plot_width*756/1024));
fprintf(fp,"%d %d %d\n",plot_height,plot_width,8);
fprintf(fp,"[%d 0 0 %d 0 %d]\n",plot_height,-plot_width,
                                plot_width);
fprintf(fp,"{currentfile redstr readhexstring pop}\n");
fprintf(fp,"{currentfile grestr readhexstring pop}\n");
fprintf(fp,"{currentfile blustr readhexstring pop}\n");
fprintf(fp,"true 3\n");
fprintf(fp,"colorimage\n");

  for ( i = plot_width-1 ; i > -1 ; i-- )
  {
    if (xinterrupt_()==1) 
    {
    fclose(fp);
    return;
    }
  XFillRectangle(mydisplay,title_window,title_gc,
                 98+i-1,0,2,20);
    for (j=0;j<plot_height;j++) 
    print_color[j].pixel=XGetPixel(image,i,j);
  XQueryColors (mydisplay,cmap,print_color,plot_height);
    for (j=0;j<plot_height;j++)
    grey[j]=print_color[j].red/256;
  pp(grey,plot_height);
    for (j=0;j<plot_height;j++)
    grey[j]=print_color[j].green/256;
  pp(grey,plot_height);
    for (j=0;j<plot_height;j++)
    grey[j]=print_color[j].blue/256;
  pp(grey,plot_height);
  }

fprintf(fp,"\nstroke\n");
fprintf(fp,"showpage\n");
fprintf(fp,"grestore\n");
fclose(fp);

XDestroyImage (image);
XClearWindow(mydisplay,title_window);

}

/* -------------------------------------------------- */
pp(g,n)

int     g[],n;

{

char     line[80];
int      line_count,i,offset;
         
line_count=0;
offset=print_count;

while (True)

{
                               
line_count++;
print_count++;

  if (print_count == 40) 
  {
  
    for (i=line_count-40+offset; i<line_count; i++)
    fprintf(fp,"%c%c",he[g[i]*2],he[g[i]*2+1]);
  fprintf(fp,"\n");
  offset=0;
  print_count=0;
  }

  if (line_count == n) 
  {
    for (i=line_count-print_count; i<line_count; i++)
    fprintf(fp,"%c%c",he[g[i]*2],he[g[i]*2+1]);                  
  return;
  }

}

}

/* -------------------------------------------------- */
xoutput_(string,nstring,dummy)

/* subroutine to print a message in the "output" window
   at the bottom left of the plot window

   arguments: string: string to be printed
              nstring: length of string
*/

char          *string;
int           *nstring, *dummy;

{      

int            i;

for (i=0;i< *nstring;i++) output_message[i]=string[i];
len_output_message= *nstring;

XClearWindow(mydisplay,output_window);
XDrawString(mydisplay,output_window,output_gc,
            2,12,string,*nstring);

}
                                                        
/* -------------------------------------------------- */
xinput_(string,nstring,dummy)

/* subroutine that prompts the user for a string to be typed
   in the input window at the bottom right of the plot window

   arguments: string: string typed in by the user
              nstring: length of string
*/

char          *string;
int           *nstring, *dummy;

{

XEvent        event;
int           count;
char          ch[20];
KeySym        key;
int           i;
int           offset;
int           status;
Time          time_focus;

*nstring=0;
offset=4;
XWarpPointer(mydisplay,None,input_window,
             0,0,0,0,
             offset,8);


while (True)
  {
  XNextEvent (mydisplay, &event);

  switch (event.type)

    {

    case KeyPress:

      if (event.xkey.window == input_window)
      count=XLookupString ( &event,ch,20, &key,0);

      {
        switch (key)

        {

        case XK_Return:
          XClearWindow(mydisplay,input_window);
          string[ *nstring]=' ';
          return;
          break;
                               
        case XK_Delete:
          *nstring= *nstring-1;
          if ( *nstring < 0 ) *nstring=0;
          break;
        
        case XK_BackSpace:
          *nstring= *nstring-1;
          if ( *nstring < 0 ) *nstring=0;
          break;
        
        default:                 
          for (i=0;i<count;i++) string[ *nstring+i] = ch[i];
          *nstring= *nstring+count;
          break;

        }

      }

      break;

    case MappingNotify:
      XRefreshKeyboardMapping( &event);
      break;

    case Expose:

      if (event.xexpose.window == plot_window)
      XCopyArea(mydisplay,save_pixmap,
                plot_window,plot_gc,
                event.xexpose.x,event.xexpose.y,
                event.xexpose.width,event.xexpose.height,
                event.xexpose.x,event.xexpose.y);
      XClearWindow(mydisplay,output_window);
      XDrawString(mydisplay,output_window,output_gc,
                  2,12,output_message,len_output_message);
 
      break;

    }
                                             
  string[ *nstring]='\0';
  XClearWindow(mydisplay,input_window);
  offset=4+XTextWidth(mediumfont_struct,string, *nstring);
  if (event.xkey.window == input_window)
    XWarpPointer(mydisplay,None,input_window,
               0,0,0,0,
               offset,8);
  XDrawString(mydisplay,input_window,input_gc,
              2,12,string,*nstring);
                          
  }

}
/* -------------------------------------------------- */
xcircle_(xc,yc,radius)

/* to draw a circle 
   arguments: xc and yc are the coordinates of the centre of the circle
              radius is the radius
*/


float     *xc, *yc, *radius;

{

int          x,y;
unsigned int width,height;
int          angle1,angle2;
float        xp,yp,radx,rady;

if (PostScript == 'y') pcircle_ (xc,yc,radius);

  switch (scale.mapping) {

  case ('w'):
    xp = *xc; yp = *yc; radx= *radius; rady= *radius;
    break;
 
  case ('u'):
    xp=( *xc-scale.xmin)/(scale.xmax-scale.xmin)
      *(scale.xpmax-scale.xpmin)+scale.xpmin;
    yp=( *yc-scale.ymin)/(scale.ymax-scale.ymin)
      *(scale.ypmax-scale.ypmin)+scale.ypmin;
    radx= *radius/(scale.xmax-scale.xmin)*(scale.xpmax-scale.xpmin);
    rady= *radius/(scale.ymax-scale.ymin)*(scale.ypmax-scale.ypmin);
    break;

  }

if (radx < 0) radx= -radx;
if (rady < 0) rady= -rady;
x=xp-radx; y=yp-rady;
width=radx*2; height=rady*2;
angle1=0;angle2=360*64;

XDrawArc (mydisplay,plot_window,plot_gc,
          x,y,width,height,angle1,angle2);
                                          
if (save_mode == 's')
XDrawArc (mydisplay,save_pixmap,plot_gc,
          x,y,width,height,angle1,angle2);

} 

/* -------------------------------------------------- */
xcirclef_(xc,yc,radius)

/* same as xcircle but produces a filled circle */


float     *xc, *yc, *radius;

{

int          x,y;
unsigned int width,height;
int          angle1,angle2;
float        xp,yp,radx,rady;

if (PostScript == 'y') pcirclef_ (xc,yc,radius);

  switch (scale.mapping) {

  case ('w'):
    xp = *xc; yp = *yc; radx= *radius; rady= *radius;
    break;
 
  case ('u'):
    xp=( *xc-scale.xmin)/(scale.xmax-scale.xmin)
      *(scale.xpmax-scale.xpmin)+scale.xpmin;
    yp=( *yc-scale.ymin)/(scale.ymax-scale.ymin)
      *(scale.ypmax-scale.ypmin)+scale.ypmin;
    radx= *radius/(scale.xmax-scale.xmin)*(scale.xpmax-scale.xpmin);
    rady= *radius/(scale.ymax-scale.ymin)*(scale.ypmax-scale.ypmin);
    break;

  }

if (radx < 0) radx= -radx;
if (rady < 0) rady= -rady;
x=xp-radx; y=yp-rady;
width=radx*2; height=rady*2;
angle1=0;angle2=360*64;

XFillArc (mydisplay,plot_window,plot_gc,
          x,y,width,height,angle1,angle2);
                                          
if (save_mode == 's')
XFillArc (mydisplay,save_pixmap,plot_gc,
          x,y,width,height,angle1,angle2);

} 

/* -------------------------------------------------- */
xlogo_(x,y)

/* produces an ANU logo
   arguments: x and y are the coordinates of the upper left corner
              of the logo
*/

int     *x, *y;

{

GC                     logo_gc;

logo_gc=XCreateGC(mydisplay,logo_pixmap,0,0);
XSetForeground (mydisplay,logo_gc,myforeground);
XSetBackground (mydisplay,logo_gc,mybackground);
                                                   
XCopyArea (mydisplay,logo_pixmap,plot_window,logo_gc,
           0,0,100,120, *x, *y);

if (save_mode == 's')
XCopyArea (mydisplay,logo_pixmap,save_pixmap,logo_gc,
           0,0,100,120, *x, *y);

}


/* ------------------------------------------------ */
xraw_(fnme,nlen,nfnme)

/* stores the content of the plot window in a raw format image
   note that this format is internal to xplot and should only
   be used to be handled with the unraw command
   arguments: fnme: name of the raw format file
              nfnme: length of file name
*/

char     *fnme;
int      *nfnme;
int      *nlen;

{
                               
XImage    *image;
int        i,j;
char       file[256];
XColor     current_color[256];
char       col[256];
static int def[]={30000,65535};

xpopup_();   

image=XGetImage (mydisplay,save_pixmap,0,0,
                 plot_width,plot_height,AllPlanes,ZPixmap);

for ( i=0 ; i< *nlen ; i++) file[i]=fnme[i];
file[ *nlen]='.';file[ *nlen+1]='r';file[ *nlen+2]='a';
file[ *nlen+3]='w';file[ *nlen+4]='\0';

fp=fopen(file,"w");

  for ( j=0 ; j<plot_height ; j++ ) 
  {
    if (xinterrupt_()==1) 
    {
    fclose(fp);
    return;
    }
  XFillRectangle(mydisplay,title_window,title_gc,
                 98+(j*plot_width)/plot_height-1,0,2,20);
    for ( i=0 ; i<plot_width ; i++ ) 
    putc(XGetPixel(image,i,j),fp);
  }

fclose(fp);

XDestroyImage (image);
XClearWindow(mydisplay,title_window);

for ( i=0 ; i< *nlen ; i++) file[i]=fnme[i];
file[ *nlen]='.';file[ *nlen+1]='p';file[ *nlen+2]='a';
file[ *nlen+3]='l';file[ *nlen+4]='\0';

fp=fopen(file,"w");
                                        
for ( i=0 ; i<128 ; i++) for ( j=0 ; j<2 ; j++) 
     col[i*2+j]=def[j]/256;
XQueryColors (mydisplay,cmap,color,ncolors);

  for (i=0 ; i<ncolors ; i++) col[color[i].pixel]=color[i].red/256;
  for (i=0 ; i<256 ; i++) putc(col[i],fp);

  for (i=0 ; i<ncolors ; i++) col[color[i].pixel]=color[i].green/256;
  for (i=0 ; i<256 ; i++) putc(col[i],fp);

  for (i=0 ; i<ncolors ; i++) col[color[i].pixel]=color[i].blue/256;
  for (i=0 ; i<256 ; i++) putc(col[i],fp);

fclose(fp);

}

/*------------------------------------------------------*/
xunraw_(fnme,nlen,nfnme)

/* brings an image saved with the raw command back to the
   plot window
   arguments: fnme: name of the raw format file
              nfnme: length of file name
*/

char   *fnme;
int    *nfnme;
int    *nlen;

{

char           file[256];
XImage        *image;
int            i,j;
char           col[256];

for ( i=0 ; i< *nlen ; i++) file[i]=fnme[i];
file[ *nlen]='.';file[ *nlen+1]='p';file[ *nlen+2]='a';
file[ *nlen+3]='l';file[ *nlen+4]='\0';

fp=fopen(file,"r");

if (fp == NULL) { *nlen=0; return;}
                                        
for (i=0 ; i<256 ; i++) col[i]=getc(fp);
for (i=0 ; i<ncolors ; i++) color[i].red=col[pixels[i]]*256;

for (i=0 ; i<256 ; i++) col[i]=getc(fp);
for (i=0 ; i<ncolors ; i++) color[i].green=col[pixels[i]]*256;

for (i=0 ; i<256 ; i++) col[i]=getc(fp);
for (i=0 ; i<ncolors ; i++) color[i].blue=col[pixels[i]]*256;

for (i=0 ; i<ncolors ; i++) 
  {
  color[i].pixel=pixels[i];
  color[i].flags= DoRed | DoGreen | DoBlue ;
  }

fclose(fp);

file[ *nlen]='.';file[ *nlen+1]='r';file[ *nlen+2]='a';
file[ *nlen+3]='w';file[ *nlen+4]='\0';

fp=fopen(file,"r");

if (fp==NULL) { *nlen=0; return;}

image=XGetImage (mydisplay,save_pixmap,0,0,
                 plot_width,plot_height,AllPlanes,ZPixmap);

  for ( j=0 ; j<plot_height ; j++ ) 
    for ( i=0 ; i<plot_width ; i++ ) 
    XPutPixel(image,i,j, getc(fp) );

XStoreColors (mydisplay,cmap,color,ncolors);

XPutImage(mydisplay,plot_window,plot_gc,image,
          0,0,0,0,plot_width,plot_height);

XDestroyImage (image);

fclose(fp);

}

/* ------------------------------------------------ */
xcmap_view_()

/* to view the color map */

{

Window      color_window;
XSizeHints  color_hint;
GC          color_gc;
int         i,j;
int         xp,yp1,yp2;
short       x,y;
int         count[256],count_max;
XImage     *image;
int         wait;
XEvent      color_event;

color_hint.x= 100;
color_hint.y= 50;
color_hint.width= 256;
color_hint.height= 300;
color_hint.flags=PPosition|PSize;

color_window=XCreateSimpleWindow(mydisplay,
  DefaultRootWindow (mydisplay),
  color_hint.x,color_hint.y,color_hint.width,color_hint.height,
  2,myforeground,mybackground);

color_gc=XCreateGC(mydisplay,color_window,0,0);

XSelectInput(mydisplay,color_window,
  ButtonPressMask|ExposureMask|PointerMotionMask);

XMapRaised(mydisplay,color_window);

  for (i=0;i<ncolors;i++)
  {
  XSetForeground (mydisplay,color_gc,color[i].pixel);
  xp = color[i].pixel ;
  yp1 = 0 ;
  yp2 = 99;
  XDrawLine (mydisplay,color_window,color_gc,
             xp,yp1,xp,yp2);
  }

image=XGetImage (mydisplay,save_pixmap,0,0,
                 plot_width,plot_height,AllPlanes,ZPixmap);
     
  for (i=0;i<256;i++)
  count[i]=0;
  
  for ( i = 0 ; i < plot_width ; i++)
  {
    for ( j = 0 ; j < plot_height ; j++) 
    count[XGetPixel(image,i,j)]++;
  }

XDestroyImage (image);

count_max=0;
  for (i=0;i<256;i++)
  {
  if (count[i] > count_max) count_max=count[i];
  }
   
XSetForeground (mydisplay,color_gc,myforeground);
  for (i=0;i<ncolors;i++)
  {
  xp = color[i].pixel;                         
  yp1=299;
  yp2=299-(count[i]*198)/count_max;
  XDrawLine (mydisplay,color_window,color_gc,
            xp,yp1,xp,yp2);
  }
xp=100;
yp1=0;
yp2=255;
XDrawLine (mydisplay,color_window,color_gc,yp1,xp,yp2,xp);

wait=0; 
while (wait==0)
  {

  XNextEvent (mydisplay,&color_event);

  switch (color_event.type){

  case Expose:
    if (color_event.xexpose.window == plot_window)
    XCopyArea(mydisplay,save_pixmap,
              plot_window,plot_gc,
              color_event.xexpose.x,color_event.xexpose.y,
              color_event.xexpose.width,color_event.xexpose.height,
              color_event.xexpose.x,color_event.xexpose.y);
    if (color_event.xexpose.window == output_window)
    {
    XClearWindow(mydisplay,output_window);
    XDrawString(mydisplay,output_window,output_gc,
                2,20,output_message,len_output_message);
    }
    if (color_event.xexpose.window == color_window)
    {
      for (i=0;i<ncolors;i++)
      {
      XSetForeground (mydisplay,color_gc,color[i].pixel);
      xp = color[i].pixel;
      yp1=0;
      yp2=99;
      XDrawLine (mydisplay,color_window,color_gc,
                 xp,yp1,xp,yp2);
      }  

    XSetForeground (mydisplay,color_gc,myforeground);
      for (i=0;i<ncolors;i++)
      {
      xp = color[i].pixel;                 
      yp1=299;
      yp2=299-(count[i]*198)/count_max;
      XDrawLine (mydisplay,color_window,color_gc,
                xp,yp1,xp,yp2);
      }
    xp=100;
    yp1=0;
    yp2=255;
    XDrawLine (mydisplay,color_window,color_gc,yp1,xp,yp2,xp);
    }
    break;

  case MappingNotify:
    XRefreshKeyboardMapping( &color_event);
    break;

  case ButtonPress:
    x = color_event.xbutton.x; y = color_event.xbutton.y;
    wait=1;
    XFreeGC(mydisplay,color_gc);
    XDestroyWindow(mydisplay,color_window);
    break;

    }

  }

}

/* ------------------------------------------------ */
xarea_(xo1,yo1,xo2,yo2) 

/* prompts the user to defined a rectangular area by
   selecting two points in the plot window
   arguments: xo1,yo1: returned coordinates of first point selected
              xo2,yo2: returned coordinates of second point selected
*/

float  *xo1, *yo1, *xo2, *yo2;

{

GC      mask_gc;
XEvent  event;
short   first;
int     x1,y1,width,height;

first=0;

mask_gc=XCreateGC (mydisplay,plot_window,0,0);
XSetFunction (mydisplay,mask_gc,GXinvert);
XSetPlaneMask (mydisplay,mask_gc,AllPlanes);

while (True)
  {
  XNextEvent (mydisplay, &event);

  switch (event.type)

    {

      case Expose:
      if (event.xexpose.window == plot_window)
      XCopyArea(mydisplay,save_pixmap,
                plot_window,plot_gc,
                event.xexpose.x,event.xexpose.y,
                event.xexpose.width,event.xexpose.height,
                event.xexpose.x,event.xexpose.y);
      if (event.xexpose.window == output_window)
      {
      XClearWindow(mydisplay,output_window);
      XDrawString(mydisplay,output_window,output_gc,
                  2,20,output_message,len_output_message);
      }
      break;

    case MappingNotify:
      XRefreshKeyboardMapping( &event);
      break;

    case MotionNotify:
      if (event.xmotion.window == plot_window & first == 1)
      {                                        
      XDrawRectangle(mydisplay,plot_window,mask_gc,x1,y1,width,height);
      width = event.xmotion.x-x1; height = event.xmotion.y-y1;
      XDrawRectangle(mydisplay,plot_window,mask_gc,x1,y1,width,height);
      }
      break;

    case ButtonPress:
      if (first == 1) 
      {
        switch (scale.mapping) {
        case ('w'):
          *xo1 = x1; *yo1 = y1;
          *xo2 = x1+width; *yo2 = y1+height;
          break;
        case ('u'):
          *xo1=((x1-scale.xpmin)
               *(scale.xmax-scale.xmin))
               /(scale.xpmax-scale.xpmin)+scale.xmin;
          *yo1=((y1-scale.ypmin)
               *(scale.ymax-scale.ymin))
               /(scale.ypmax-scale.ypmin)+scale.ymin;
          *xo2=((x1+width-scale.xpmin)
               *(scale.xmax-scale.xmin))
               /(scale.xpmax-scale.xpmin)+scale.xmin;
          *yo2=((y1+height-scale.ypmin)
               *(scale.ymax-scale.ymin))
               /(scale.ypmax-scale.ypmin)+scale.ymin;
          break;
          }                     
      return;
      }

      if (first == 0) 
      {
      x1=event.xmotion.x;
      y1=event.xmotion.y;
      first=1;
      width=1;
      height=1;
      XDrawRectangle(mydisplay,plot_window,mask_gc,x1,y1,width,height);
      }
      break;

    }

  }

}

/* ------------------------------------------------ */
xsection_(xo1,yo1,xo2,yo2) 

 /* prompts the user to defined a line by
   selecting two points in the plot window
   arguments: xo1,yo1: returned coordinates of first point selected
              xo2,yo2: returned coordinates of second point selected
*/
                     
float  *xo1, *yo1, *xo2, *yo2;

{

GC      mask_gc;
XEvent  event;
short   first;
int     x1,y1,width,height;

first=0;

mask_gc=XCreateGC (mydisplay,plot_window,0,0);
XSetFunction (mydisplay,mask_gc,GXinvert);
XSetPlaneMask (mydisplay,mask_gc,AllPlanes);

while (True)
  {
  XNextEvent (mydisplay, &event);

  switch (event.type)

    {

      case Expose:
      if (event.xexpose.window == plot_window)
      XCopyArea(mydisplay,save_pixmap,
                plot_window,plot_gc,
                event.xexpose.x,event.xexpose.y,
                event.xexpose.width,event.xexpose.height,
                event.xexpose.x,event.xexpose.y);
      if (event.xexpose.window == output_window)
      {
      XClearWindow(mydisplay,output_window);
      XDrawString(mydisplay,output_window,output_gc,
                  2,20,output_message,len_output_message);
      }
      break;

    case MappingNotify:
      XRefreshKeyboardMapping( &event);
      break;

    case MotionNotify:
      if (event.xmotion.window == plot_window & first == 1)
      {                                        
      XDrawLine(mydisplay,plot_window,mask_gc,x1,y1,x1+width,y1+height);
      width = event.xmotion.x-x1; height = event.xmotion.y-y1;
      XDrawLine(mydisplay,plot_window,mask_gc,x1,y1,x1+width,y1+height);
      }
      break;

    case ButtonPress:
      if (first == 1) 
      {
        switch (scale.mapping) {
        case ('w'):
          *xo1 = x1; *yo1 = y1;
          *xo2 = x1+width; *yo2 = y1+height;
          break;
        case ('u'):
          *xo1=((x1-scale.xpmin)
               *(scale.xmax-scale.xmin))
               /(scale.xpmax-scale.xpmin)+scale.xmin;
          *yo1=((y1-scale.ypmin)
               *(scale.ymax-scale.ymin))
               /(scale.ypmax-scale.ypmin)+scale.ymin;
          *xo2=((x1+width-scale.xpmin)
               *(scale.xmax-scale.xmin))
               /(scale.xpmax-scale.xpmin)+scale.xmin;
          *yo2=((y1+height-scale.ypmin)
               *(scale.ymax-scale.ymin))
               /(scale.ypmax-scale.ypmin)+scale.ymin;
          break;
          }                     
      return;
      }

      if (first == 0) 
      {
      x1=event.xmotion.x;
      y1=event.xmotion.y;
      first=1;
      width=1;
      height=1;
      XDrawLine(mydisplay,plot_window,mask_gc,x1,y1,x1+width,y1+height);
      }
      break;

    }

  }

}

/* ------------------------------------------------ */
xpolygon_(xpol,ypol,npol)

 /* prompts the user to defined a polygon by
   selecting several points in the plot window
   arguments: xpol,ypol: returned coordinates of corners of the
                         polygon
              npol: number of corners to the polygon
   note that the user should use the left mouse button to select the
   corners and the right mouse button to finish
*/
                     
int    *npol;
float  *xpol, *ypol;

{

GC      mask_gc;
XEvent  event;
XPoint  points[1000];
int     npoints,dx,dy,i;
int     xm,ym;
               
npoints=0;
xm=0;ym=0;

XDefineCursor(mydisplay,plot_window,arrowcursor);
mask_gc=XCreateGC (mydisplay,plot_window,0,0);
XSetFunction (mydisplay,mask_gc,GXinvert);
XSetPlaneMask (mydisplay,mask_gc,AllPlanes);

XDrawLine(mydisplay,plot_window,mask_gc,
          xm,0,xm,plot_height);
XDrawLine(mydisplay,plot_window,mask_gc,
          0,ym,plot_width,ym);

while (True)
  {
  XNextEvent (mydisplay, &event);

  switch (event.type)

    {

      case Expose:
      if (event.xexpose.window == plot_window)
      XCopyArea(mydisplay,save_pixmap,
                plot_window,plot_gc,
                event.xexpose.x,event.xexpose.y,
                event.xexpose.width,event.xexpose.height,
                event.xexpose.x,event.xexpose.y);
      if (event.xexpose.window == output_window)
      {
      XClearWindow(mydisplay,output_window);
      XDrawString(mydisplay,output_window,output_gc,
                  2,20,output_message,len_output_message);
      }
      break;

    case MappingNotify:
      XRefreshKeyboardMapping( &event);
      break;

    case MotionNotify:
      XDrawLine(mydisplay,plot_window,mask_gc,
                xm,0,xm,plot_height);
      XDrawLine(mydisplay,plot_window,mask_gc,
                0,ym,plot_width,ym);
      xm=event.xmotion.x; ym=event.xmotion.y;
      XDrawLine(mydisplay,plot_window,mask_gc,
                xm,0,xm,plot_height);
      XDrawLine(mydisplay,plot_window,mask_gc,
                0,ym,plot_width,ym);
      if (event.xmotion.window == plot_window & npoints != 0)
      {                                        
      XDrawLines(mydisplay,plot_window,mask_gc,
                 points,npoints+1,CoordModeOrigin);
      points[npoints].x=event.xmotion.x;
      points[npoints].y=event.xmotion.y;
      XDrawLines(mydisplay,plot_window,mask_gc,
                 points,npoints+1,CoordModeOrigin);
      }
      break;

    case ButtonPress:
      XDrawLines(mydisplay,plot_window,mask_gc,
                 points,npoints+1,CoordModeOrigin);
      points[npoints].x=event.xmotion.x;
      points[npoints].y=event.xmotion.y;
      npoints++;
      points[npoints].x=event.xmotion.x;
      points[npoints].y=event.xmotion.y;
      XDrawLines(mydisplay,plot_window,mask_gc,
                 points,npoints+1,CoordModeOrigin);
      dx=points[npoints-1].x-points[0].x;
      dy=points[npoints-1].y-points[0].y;
      if ( ((dx*dx+dy*dy) < 25  && npoints != 1) || npoints == 1000
          || event.xbutton.button == 3) 
      {
      *npol = npoints-1;
        switch (scale.mapping)
        {
        case ('w'):
          for (i=0 ; i<npoints-1 ; i++)
          {
          xpol[i]=points[i].x;
          ypol[i]=points[i].y;
          }
          break;
        case ('u'):
          for (i=0 ; i<npoints-1 ; i++)
          {
          xpol[i]=((points[i].x-scale.xpmin)
                  *(scale.xmax-scale.xmin))
                  /(scale.xpmax-scale.xpmin)+scale.xmin;
          ypol[i]=((points[i].y-scale.ypmin)
                  *(scale.ymax-scale.ymin))
                  /(scale.ypmax-scale.ypmin)+scale.ymin;
          }
          break;
        }                     
      XDrawLine(mydisplay,plot_window,mask_gc,
                xm,0,xm,plot_height);
      XDrawLine(mydisplay,plot_window,mask_gc,
                0,ym,plot_width,ym);
      XDefineCursor(mydisplay,plot_window,mycursor);
      return;
      }
      break;

    }

  }

}

/* ------------------------------------------------ */
xshapepol_(xpol,ypol,npol)

/* to reshape a polygon by either moving the corners or moving
   the whole polygon
   arguments: xpol, ypol are the coordinates of the corners of the
                         polygon
              npol is the number of corners to the polygon
*/
                      
int    *npol;
float  *xpol, *ypol;

{

GC      mask_gc;
XEvent  event;
XPoint  points[1000];
int     npoints,dx,dy,i;
int     xm,ym,corner;
int     xc,yc;
               
npoints= *npol;
switch (scale.mapping)
  {
  case ('w'):
    for (i=0 ; i<npoints ; i++)
    {
    points[i].x=xpol[i];
    points[i].y=ypol[i];
    }
    break;
  case ('u'):
    for (i=0 ; i<npoints ; i++)
    {
    points[i].x=((xpol[i]-scale.xmin)
               *(scale.xpmax-scale.xpmin))
               /(scale.xmax-scale.xmin)+scale.xpmin;
    points[i].y=((ypol[i]-scale.ymin)
               *(scale.ypmax-scale.ypmin))
               /(scale.ymax-scale.ymin)+scale.ypmin;
    }
  points[npoints].x=points[0].x;  
  points[npoints].y=points[0].y;  
  break;
  }                     
                         
xc=0;yc=0;
for (i=0;i<npoints;i++) {xc=xc+points[i].x;yc=yc+points[i].y;}
xc=xc/npoints;yc=yc/npoints;

xm=0;ym=0;corner= -1;

XDefineCursor(mydisplay,plot_window,arrowcursor);
mask_gc=XCreateGC (mydisplay,plot_window,0,0);
XSetFunction (mydisplay,mask_gc,GXinvert);
XSetPlaneMask (mydisplay,mask_gc,AllPlanes);

XDrawLine(mydisplay,plot_window,mask_gc,
          xm,0,xm,plot_height);
XDrawLine(mydisplay,plot_window,mask_gc,
          0,ym,plot_width,ym);

XDrawLine(mydisplay,plot_window,mask_gc,
          xc-5,yc,xc+5,yc);
XDrawLine(mydisplay,plot_window,mask_gc,
          xc,yc-5,xc,yc+5);
XDrawLines(mydisplay,plot_window,mask_gc,
           points,npoints+1,CoordModeOrigin);

while (True)
  {
  XNextEvent (mydisplay, &event);

  switch (event.type)

    {

    case Expose:
      if (event.xexpose.window == plot_window)
      XCopyArea(mydisplay,save_pixmap,
                plot_window,plot_gc,
                event.xexpose.x,event.xexpose.y,
                event.xexpose.width,event.xexpose.height,
                event.xexpose.x,event.xexpose.y);
      if (event.xexpose.window == output_window)
      {
      XClearWindow(mydisplay,output_window);
      XDrawString(mydisplay,output_window,output_gc,
                  2,20,output_message,len_output_message);
      }
      break;

    case MappingNotify:
      XRefreshKeyboardMapping( &event);
      break;

    case MotionNotify:
      XDrawLine(mydisplay,plot_window,mask_gc,
                xm,0,xm,plot_height);
      XDrawLine(mydisplay,plot_window,mask_gc,
                0,ym,plot_width,ym);
      xm=event.xmotion.x; ym=event.xmotion.y;
      XDrawLine(mydisplay,plot_window,mask_gc,
                xm,0,xm,plot_height);
      XDrawLine(mydisplay,plot_window,mask_gc,
                0,ym,plot_width,ym);
      if (event.xmotion.window == plot_window & corner >= 0)
      {                                        
      XDrawLine(mydisplay,plot_window,mask_gc,
                xc-5,yc,xc+5,yc);
      XDrawLine(mydisplay,plot_window,mask_gc,
                xc,yc-5,xc,yc+5);
      XDrawLines(mydisplay,plot_window,mask_gc,
                 points,npoints+1,CoordModeOrigin);
      points[corner].x=event.xmotion.x;
      points[corner].y=event.xmotion.y;
      points[npoints].x=points[0].x;
      points[npoints].y=points[0].y;
      xc=0;yc=0;
      for (i=0;i<npoints;i++) {xc=xc+points[i].x;yc=yc+points[i].y;}
      xc=xc/npoints;yc=yc/npoints;
      XDrawLine(mydisplay,plot_window,mask_gc,
                xc-5,yc,xc+5,yc);
      XDrawLine(mydisplay,plot_window,mask_gc,
                xc,yc-5,xc,yc+5);
      XDrawLines(mydisplay,plot_window,mask_gc,
                 points,npoints+1,CoordModeOrigin);
      }
      if (event.xmotion.window == plot_window & corner == -2)
      {                                        
      XDrawLine(mydisplay,plot_window,mask_gc,
                xc-5,yc,xc+5,yc);
      XDrawLine(mydisplay,plot_window,mask_gc,
                xc,yc-5,xc,yc+5);
      XDrawLines(mydisplay,plot_window,mask_gc,
                 points,npoints+1,CoordModeOrigin);
      dx=event.xmotion.x-xc;
      dy=event.xmotion.y-yc;
        for (i=0;i<npoints;i++)
        {
        points[i].x=points[i].x+dx;
        points[i].y=points[i].y+dy;
        }
      points[npoints].x=points[0].x;
      points[npoints].y=points[0].y;
      xc=0;yc=0;
      for (i=0;i<npoints;i++) {xc=xc+points[i].x;yc=yc+points[i].y;}
      xc=xc/npoints;yc=yc/npoints;
      XDrawLine(mydisplay,plot_window,mask_gc,
                xc-5,yc,xc+5,yc);
      XDrawLine(mydisplay,plot_window,mask_gc,
                xc,yc-5,xc,yc+5);
      XDrawLines(mydisplay,plot_window,mask_gc,
                 points,npoints+1,CoordModeOrigin);
      }
      break;

    case ButtonPress:
      switch (corner)
        {
        case -1:
        XDrawLine(mydisplay,plot_window,mask_gc,
                  xc-5,yc,xc+5,yc);
        XDrawLine(mydisplay,plot_window,mask_gc,
                  xc,yc-5,xc,yc+5);
        XDrawLines(mydisplay,plot_window,mask_gc,
           points,npoints+1,CoordModeOrigin);
          for (i=0;i<npoints;i++)
          {
          dx=points[i].x-event.xmotion.x;
          dy=points[i].y-event.xmotion.y;
          if ( (dx*dx+dy*dy) < 25  ) corner=i;
          }
        dx=xc-event.xmotion.x;
        dy=yc-event.xmotion.y;
        if ( (dx*dx+dy*dy) < 25  ) corner= -2;
        XDrawLine(mydisplay,plot_window,mask_gc,
                  xc-5,yc,xc+5,yc);
        XDrawLine(mydisplay,plot_window,mask_gc,
                  xc,yc-5,xc,yc+5);
        XDrawLines(mydisplay,plot_window,mask_gc,
                   points,npoints+1,CoordModeOrigin);
        break;
        default:
        switch (scale.mapping)
          {
          case ('w'):
            for (i=0 ; i<npoints ; i++)
            {
            xpol[i]=points[i].x;
            ypol[i]=points[i].y;
            }
            break;
          case ('u'):
            for (i=0 ; i<npoints ; i++)
            {
            xpol[i]=((points[i].x-scale.xpmin)
                    *(scale.xmax-scale.xmin))
                    /(scale.xpmax-scale.xpmin)+scale.xmin;
            ypol[i]=((points[i].y-scale.ypmin)
                    *(scale.ymax-scale.ymin))
                    /(scale.ypmax-scale.ypmin)+scale.ymin;
            }
            break;
          }                     
        XDrawLine(mydisplay,plot_window,mask_gc,
                  xm,0,xm,plot_height);
        XDrawLine(mydisplay,plot_window,mask_gc,
                  0,ym,plot_width,ym);
        XDrawLines(mydisplay,plot_window,mask_gc,
                   points,npoints+1,CoordModeOrigin);
        XDefineCursor(mydisplay,plot_window,mycursor);
        return;
        break;
        }
    break;

    }

  }

}

/* ------------------------------------------------- */
xflush_()

/* purpose: to flush the buffer
*/

{

if (PostScript == 'y') pflush_ ();

XFlush(mydisplay);
                                                                
}

/* ------------------------------------------------ */
xphoto_(fnme,nlen,nfnme)

/* dumps the content of the plot window in a file to be
   read as a "raw" format image in photoshop

   arguments: fnme: name of the file
              nlen: length of name of the file
*/

char     *fnme;
int      *nfnme;
int      *nlen;

{
                               
XImage    *image;
int        i,j;
char       file[256];
XColor     print_color[1200];
char       index;

xpopup_();   

image=XGetImage (mydisplay,save_pixmap,0,0,
                 plot_width,plot_height,AllPlanes,ZPixmap);

for ( i=0 ; i< *nlen ; i++) file[i]=fnme[i];
file[ *nlen]='\0';

fp=fopen(file,"w");
  if (fp==NULL) 
  {
  printf("Could not open file %s\n",file);
  return;
  }

  for ( j=0 ; j<plot_height ; j++ ) 
  {
    if (xinterrupt_()==1) 
    {
    fclose(fp);
    return;
    }
  XFillRectangle(mydisplay,title_window,title_gc,
                 98+(j*plot_width)/plot_height-1,0,2,20);
    for ( i=0 ; i<plot_width ; i++ ) 
    print_color[i].pixel=XGetPixel(image,i,j);
  XQueryColors (mydisplay,cmap,print_color,plot_width);
    for ( i=0 ; i<plot_width ; i++)
    {    
    index=(.35*print_color[i].red+
           .55*print_color[i].green+
           .10*print_color[i].blue)/256;
    putc(index,fp);
    }
  }

fclose(fp);

XDestroyImage (image);
XClearWindow(mydisplay,title_window);

}

/* ------------------------------------------------ */
xcolorphoto_(fnme,nlen,nfnme)

/* same as xphoto but produces a true color image */

char     *fnme;
int      *nfnme;
int      *nlen;

{
                               
XImage    *image;
int        i,j;
char       file[256];
XColor     print_color[1200];

xpopup_();   

image=XGetImage (mydisplay,save_pixmap,0,0,
                 plot_width,plot_height,AllPlanes,ZPixmap);

for ( i=0 ; i< *nlen ; i++) file[i]=fnme[i];
file[ *nlen]='\0';

fp=fopen(file,"w");
  if (fp==NULL) 
  {
  printf("Could not open file %s\n",file);
  return;
  }

  for ( j=0 ; j<plot_height ; j++ ) 
  {
    if (xinterrupt_()==1) 
    {
    fclose(fp);
    return;
    }
  XFillRectangle(mydisplay,title_window,title_gc,
                 98+(j*plot_width)/plot_height-1,0,2,20);
    for ( i=0 ; i<plot_width ; i++ ) 
    print_color[i].pixel=XGetPixel(image,i,j);
  XQueryColors (mydisplay,cmap,print_color,plot_width);
    for ( i=0 ; i<plot_width ; i++)
    {
    putc(print_color[i].red/256,fp);
    putc(print_color[i].green/256,fp);
    putc(print_color[i].blue/256,fp);
    }
  }

fclose(fp);

XDestroyImage (image);
XClearWindow(mydisplay,title_window);

}

/* ------------------------------------------------ */
int xinterrupt_()

/* to check if the interrupt button has been pressed
   the interrupt button is the small button on the top left corner
   of the plot window
   xinterrupt=0 : interrupt button has not been pressed
   xinterrupt=1 : interrupt button has been pressed
*/

{

XEvent   interrupt_event;
int      iflag;

iflag=0;
if (XCheckWindowEvent(mydisplay,
                      interrupt_window,
                      ButtonPressMask,
                      &interrupt_event)) iflag=1;
return(iflag);

}

/* ------------------------------------------------ */
xp6_(fnme,nlen,nfnme)

char     *fnme;
int      *nfnme;
int      *nlen;

{
                               
XImage    *image;
int        i,j;
char       file[256];
XColor     print_color[1200];

xpopup_();   

image=XGetImage (mydisplay,save_pixmap,0,0,
                 plot_width,plot_height,AllPlanes,ZPixmap);

for ( i=0 ; i< *nlen ; i++) file[i]=fnme[i];
file[ *nlen]='\0';

fp=fopen(file,"w");
  if (fp==NULL) 
  {
  printf("Could not open file %s\n",file);
  return;
  }

fprintf(fp,"P6 %d %d 255\n",plot_width,plot_height);

  for ( j=0 ; j<plot_height ; j++ ) 
  {
    if (xinterrupt_()==1) 
    {
    fclose(fp);
    return;
    }
  XFillRectangle(mydisplay,title_window,title_gc,
                 98+(j*plot_width)/plot_height-1,0,2,20);
    for ( i=0 ; i<plot_width ; i++ ) 
    print_color[i].pixel=XGetPixel(image,i,j);
  XQueryColors (mydisplay,cmap,print_color,plot_width);
    for ( i=0 ; i<plot_width ; i++)
    {
    putc(print_color[i].red/256,fp);
    putc(print_color[i].green/256,fp);
    putc(print_color[i].blue/256,fp);
    }
  }

fclose(fp);

XDestroyImage (image);
XClearWindow(mydisplay,title_window);

}

/* ------------------------------------------------ */
xp5_(fnme,nlen,nfnme)

char     *fnme;
int      *nfnme;
int      *nlen;

{

XImage    *image;
int        i,j;
char       file[256];
XColor     print_color[1200];

xpopup_();

image=XGetImage (mydisplay,save_pixmap,0,0,
                 plot_width,plot_height,AllPlanes,ZPixmap);

for ( i=0 ; i< *nlen ; i++) file[i]=fnme[i];
file[ *nlen]='\0';

fp=fopen(file,"w");
  if (fp==NULL)
  {
  printf("Could not open file %s\n",file);
  return;
  }

fprintf(fp,"P5 %d %d 255\n",plot_width,plot_height);

  for ( j=0 ; j<plot_height ; j++ )
  {
    if (xinterrupt_()==1)
    {
    fclose(fp);
    return;
    }
  XFillRectangle(mydisplay,title_window,title_gc,
                 98+(j*plot_width)/plot_height-1,0,2,20);
    for ( i=0 ; i<plot_width ; i++ )
    print_color[i].pixel=XGetPixel(image,i,j);
  XQueryColors (mydisplay,cmap,print_color,plot_width);
    for ( i=0 ; i<plot_width ; i++)
    {
    putc((.35*print_color[i].red+
          .55*print_color[i].green+
          .10*print_color[i].blue)/256,fp);
    }
  }

fclose(fp);

XDestroyImage (image);
XClearWindow(mydisplay,title_window);

}

/*-----------------------------------------------------*/
xopen_bitmap_(sizex,sizey)

int  *sizex, *sizey;

{

interim_plot_window=plot_window;
interim_plot_gc=plot_gc;
bitmap_width= *sizex; bitmap_height= *sizey;
plot_window=XCreatePixmap (mydisplay,master_window,
                           bitmap_width,bitmap_height, 
                           depth);
if (plot_window == 0) 
 {printf("Could not create bitmap...\n");return;}
 
plot_gc=XCreateGC(mydisplay,plot_window,0,0);
XSetForeground (mydisplay,plot_gc,0);
XFillRectangle (mydisplay,plot_window,plot_gc,0,0, *sizex, *sizey);

}

/*-----------------------------------------------------*/
xclose_bitmap_()

{

XFreePixmap (mydisplay,plot_window);
plot_window=interim_plot_window;
XFreeGC(mydisplay,plot_gc);
plot_gc=interim_plot_gc;
XDestroyImage (bitmap_image);

}

/*---------------------------------------------------*/
xpixel_bitmap_ (ix,iy,bred,bgreen,bblue)

int   *ix, *iy;
int   *bred, *bgreen, *bblue;
                        
{

long               bitmap_pixel;
int                i;
    
if (bitmap_width != 0) 
{
bitmap_image=XGetImage (mydisplay,plot_window,0,0,
                        bitmap_width,bitmap_height,
                        AllPlanes,ZPixmap);
bitmap_width=0;
}
           

bitmap_pixel=XGetPixel (bitmap_image, *ix, *iy);
  if (bitmap_pixel == 0) 
  {*bred= -1; *bgreen= -1; *bblue= -1; return;}

  for (i=0;i<ncolors;i++)
  if (bitmap_pixel == color[i].pixel)
   { *bred=color[i].red;
     *bgreen=color[i].green;
     *bblue=color[i].blue;
      return;}

}

/* -------------------------------------------------- */
xpix_(x,y,p)

int    *x, *y, *p;

{
                               
XImage         *image;
unsigned long   pixel;
int             i;

image=XGetImage (mydisplay,save_pixmap,0,0,
                 plot_width,plot_height,AllPlanes,ZPixmap);

pixel=XGetPixel(image, *x, *y);
                                    
*p = 0;
  for (i=0;i<ncolors;i++)
  {
  if (pixel == color[i].pixel) *p = i;
  }

XDestroyImage (image);
}

/* -------------------------------------------------- */
xwait_()

/* purpose: to activate the wait cursor */

{

mycursor=XCreateFontCursor(mydisplay,XC_watch);

XDefineCursor(mydisplay,plot_window,mycursor);

}

/* -------------------------------------------------- */
xsupersymbol_ (x,y,size,string,nlen,dummy)

float *x, *y;
int   *size;
char  *string;
int   *nlen;
int   *dummy;

{

Font             special_font;
XFontStruct     *special_font_struct;
short            xstring,ystring;
int              width;

if (PostScript == 'y') psupersymbol_ (x,y,size,string,nlen,dummy);

 switch ( *size) {
  case 1: special_font=XLoadFont(mydisplay,"9x15");break;
  case 2: special_font=XLoadFont(mydisplay,"9x15");break;
  case 3: special_font=XLoadFont(mydisplay,"9x15");break;
  case 4: special_font=XLoadFont(mydisplay,"9x15");break;
  case 5: special_font=XLoadFont(mydisplay,"9x15");break;
  case 6: special_font=XLoadFont(mydisplay,"9x15");break;
  case 7: special_font=XLoadFont(mydisplay,"9x15");break;
  case 8: special_font=XLoadFont(mydisplay,"9x15");break;
  case -1:  special_font=XLoadFont(mydisplay,"Symbol6");break;
  case -2: special_font=XLoadFont(mydisplay,"Symbol8");break;
  case -3: special_font=XLoadFont(mydisplay,"Symbol10");break;
  case -4: special_font=XLoadFont(mydisplay,"Symbol12");break;
  case -5: special_font=XLoadFont(mydisplay,"Symbol14");break;
  case -6: special_font=XLoadFont(mydisplay,"Symbol18");break;
  case -7: special_font=XLoadFont(mydisplay,"Symbol24");break;
  case -8: special_font=XLoadFont(mydisplay,"Symbol48");break;}

XSetFont(mydisplay,plot_gc,special_font);

/* find position of string depending on mapping */

switch (scale.mapping) {

  case ('w'):
    xstring = *x; ystring = *y;
    break;
 
  case ('u'):
    xstring=( *x-scale.xmin)/(scale.xmax-scale.xmin)
         *(scale.xpmax-scale.xpmin)+scale.xpmin;
    ystring=( *y-scale.ymin)/(scale.ymax-scale.ymin)
         *(scale.ypmax-scale.ypmin)+scale.ypmin;
    break;

  }

/* draw the string */

XDrawString(mydisplay,plot_window,plot_gc,
            xstring,ystring,string,*nlen);

if (save_mode == 's') 
XDrawString(mydisplay,save_pixmap,plot_gc,
            xstring,ystring,string,*nlen);

/* modifies x to be at the end of the string */

special_font_struct=XQueryFont (mydisplay,special_font);
width=XTextWidth (special_font_struct,string, *nlen);

switch (scale.mapping) {

  case ('w'):
    *x= *x+width;
    break;
 
  case ('u'):
    *x= *x + (width*(scale.xmax-scale.xmin))/(scale.xpmax-scale.xpmin);

    break;

  }

XFreeFont (mydisplay,special_font_struct);

}
/* -------------------------------------------------- */
xname_(fnme,nlen,nfnme)

/* to give a name to the window and icon */

char     *fnme;
int      *nfnme;
int      *nlen;

{

char    name[256];
int     i;

for ( i=0 ; i< *nlen ; i++) name[i]=fnme[i];
name[ *nlen]='\0';

XStoreName (mydisplay,master_window,name);
XSetIconName (mydisplay,master_window,name);

}
/* ------------------------------------------------ */
xsuperoutput_(message,n,nlen)

/* to write n messages (of length 128) to a special (large) output window */

char    *message;
int     *n;
int     *nlen;

{

Window      super_window;
XSizeHints  super_hint;
GC          super_gc;
int         wait;
int         i,j;
XEvent      super_event;
char        super_string[128];
Cursor      super_cursor;

super_hint.x= 100;
super_hint.y= 50;
super_hint.width= 900;
super_hint.height= (*n+1)*20;
super_hint.flags=PPosition|PSize;

super_window=XCreateSimpleWindow(mydisplay,
  DefaultRootWindow (mydisplay),
  super_hint.x,super_hint.y,super_hint.width,super_hint.height,
  2,myforeground,mybackground);

super_gc=XCreateGC(mydisplay,super_window,0,0);

XSelectInput(mydisplay,super_window,
  ButtonPressMask|ExposureMask|PointerMotionMask);

XMapRaised(mydisplay,super_window);

super_cursor=XCreateFontCursor(mydisplay,XC_gumby);
XDefineCursor(mydisplay,super_window,super_cursor);
XSetForeground(mydisplay,super_gc,myforeground);
XStoreName (mydisplay,super_window,"XOutput\n");
XSetIconName (mydisplay,super_window,"XOutput\n");

XSetFont(mydisplay,super_gc,mediumfont);

for (i=0;i< *n;i++)
    {
    for (j=0;j<128;j++)
         super_string[j]= message[i*128+j];
         XDrawString(mydisplay,super_window,super_gc,2,20*i+18,
         super_string,128);
    }    

wait=0; 
while (wait==0)
  {

  XNextEvent (mydisplay,&super_event);

  switch (super_event.type){

  case Expose:
    if (super_event.xexpose.window == plot_window)
    XCopyArea(mydisplay,save_pixmap,
              plot_window,plot_gc,
              super_event.xexpose.x,super_event.xexpose.y,
              super_event.xexpose.width,super_event.xexpose.height,
              super_event.xexpose.x,super_event.xexpose.y);
    if (super_event.xexpose.window == output_window)
    {
    XClearWindow(mydisplay,output_window);
    XDrawString(mydisplay,output_window,output_gc,
                2,20,output_message,len_output_message);
    }
    if (super_event.xexpose.window == super_window)
    {
     for (i=0;i< *n;i++)
         {
         for (j=0;j<128;j++)
              super_string[j]= message[i*128+j];
              XDrawString(mydisplay,super_window,super_gc,2,20*i+18,
              super_string,128);
         }    
    }
    break;

  case MappingNotify:
    XRefreshKeyboardMapping( &super_event);
    break;

  case ButtonPress:
    wait=1;
    XFreeGC(mydisplay,super_gc);
    XDestroyWindow(mydisplay,super_window);
    break;

    }

  }

}

/* -------------------------------------------------- */
xbell_()

/* to ring the bell */

{

XBell(mydisplay,-75);

}

/*-------------------------------------------------*/
ximage_(byte,bwidth,bheight,xoff,yoff,nbyte)

/* to draw an image on the plot window
 where  byte is the address of the first byte of data making up theimage
        bwidth is the width of the image (total)
        bheight is the height of the image (total)
        xoff is the x-offset pf where the image is to be plotted
        yoff is the y-offset of where the image has to be plotted
*/

char     *byte;
int      *bwidth, *bheight, *xoff, *yoff, *nbyte;

{

XImage     *image;
int        x,y,xdest,ydest;

image=XCreateImage (mydisplay,
                    DefaultVisual(mydisplay,myscreen),depth,ZPixmap,0,byte,
                    *bwidth, *bheight,32,0);

x=0; if ( *xoff>0) x= *xoff;
y=0; if ( *yoff>0) y= *yoff;

xdest=0; if ( *xoff<0) xdest=- *xoff;
ydest=0; if ( *yoff<0) ydest=- *yoff;

XPutImage(mydisplay,plot_window,plot_gc,image,x,y,xdest,ydest,
          *bwidth, *bheight);

}
/* -------------------------------------------------- */
xkey_(returnedkey,nstring)

/* returns an ascii character if a key has been pressed recently */

char    *returnedkey;
int     *nstring;

{

XEvent        event;
int           count;
char          ch[1];
KeySym        key;
int           i;

while (True)
  {
  XNextEvent (mydisplay, &event);

  switch (event.type)

    {

    case KeyPress:

      if (event.xkey.window == master_window)
      {
        count=XLookupString ( &event,ch,1, &key,0);
        switch (key)

        {

        case XK_Up:
          returnedkey[0]=240;return;
          break;

        case XK_Home:
          returnedkey[0]=241;return;
          break;

        case XK_Left:
          returnedkey[0]=242;return;
          break;

        case XK_End:
          returnedkey[0]=243;return;
          break;

        case XK_Down:
          returnedkey[0]=244;return;
          break;

        case XK_Next:
          returnedkey[0]=245;return;
          break;

        case XK_Right:
          returnedkey[0]=246;return;
          break;

        case XK_Prior:
          returnedkey[0]=247;return;
          break;

        default:                 
          returnedkey[0]=ch[0];return;
          break;

        }
      }

      break;

    case MappingNotify:
      XRefreshKeyboardMapping( &event);
      break;

    case Expose:

      if (event.xexpose.window == plot_window)
      XCopyArea(mydisplay,save_pixmap,
                plot_window,plot_gc,
                event.xexpose.x,event.xexpose.y,
                event.xexpose.width,event.xexpose.height,
                event.xexpose.x,event.xexpose.y);
      XClearWindow(mydisplay,output_window);
      XDrawString(mydisplay,output_window,output_gc,
                  2,12,output_message,len_output_message);
 
      break;

    }
                                             
  }

}

/*-----------------------------*/
xdot_ (xm,ym)

int     *xm, *ym;

{

GC      mask_gc; 

mask_gc=XCreateGC (mydisplay,plot_window,0,0);
XSetFunction (mydisplay,mask_gc,GXinvert);
XSetPlaneMask (mydisplay,mask_gc,AllPlanes);

XDrawLine(mydisplay,plot_window,mask_gc,
          *xm-5, *ym, *xm+5, *ym);
XDrawLine(mydisplay,plot_window,mask_gc,
          *xm, *ym-5, *xm, *ym+5);

}
