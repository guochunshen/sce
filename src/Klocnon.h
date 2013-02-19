#define max_no_t 100
#define max_ext_t 100
#define pi 3.141592653589793238462643
#define SQR(a) ((a)*(a))

#define CHECK 0

typedef struct point *Pointptr;
typedef struct point { double x, y; int number; Pointptr next; } Point;

int max_no_point;

double min(double a1, double a2);
Point xy(int i, int j);
double distance(Point x1, Point x2);
void init(Pointptr *nearest);
void insert(int i, int j, Pointptr *nearest, int interval_no, Pointptr point);
void empty_nearest(Pointptr *nearest, int ant_pkt);
double epker(double t, double h);
double edgecorrection(Point p, double h);
double lm(Point p, double *par, int no, int parametric);
double transinv_edgecorrection(Point p1, Point p2);
double transinv_edgecorrection_weed(Point p1, Point p2);
double transinv_edgecorrection_rectangle(Point p1, Point p2, double sideEW, double sideNS);
double transinv_edgecorrection_spruce(Point p1, Point p2);
void gtrans(int *ant_pkt_dat, double *pointsx, double *pointsy, double *par, double *max_t, double *g_estimatedat, double *bandwidth, int *correctype, double *sideEW, double *sideNS, int *parametric, int *adaptive, int *kerneltype);
double edge_zero_correction(double t, double bw);
double neighbours2(int i, int interval_no, Pointptr *nearest, Point *points, double lm_estimate[], int correctype, double sideEW, double sideNS);
void Ktrans(int *ant_pkt_dat, double *pointsx, double *pointsy, double *par, double *max_t, double *K_estimatedat, int *correctype, double *sideEW, double *sideNS, int *parametric);
double edge(Point x, Point y, int correctype, double sideEW, double sideNS);
double neighbours2_cross(int i, int interval_no, Pointptr *nearest, Point *pointsi, double lm_estimatei[], double lm_estimatej[], int correctype, double sideEW, double sideNS);
void Ktrans_cross(int *ant_pkt_dati, double *pointsxi, double *pointsyi, int *ant_pkt_datj, double *pointsxj, double *pointsyj, double *pari, double *parj, double *max_t, double *K_estimatedat, int *correctype, double *sideEW, double *sideNS, int *parametric);
double kernel(double r, double d, double bw, int adaptive, int type);
double uniform(double t, double h);
double k1(double d1, double d2, double dist);
double k2(double d1, double d2, double dist);
double ripley_edgecorrection(Point p1, Point p2, double sideEW, double sideNS);
