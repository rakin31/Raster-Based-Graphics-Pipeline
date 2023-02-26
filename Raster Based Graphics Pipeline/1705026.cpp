#include<iostream>
#include<fstream>
#include<stack>
#include<string>
#include<math.h>
#include<time.h>
#include<cstdlib>
#include "bitmap_image.hpp"

#define pi (2*acos(0.0))
#define inf 9999999

using namespace std;


int push_count = 0;
int triangle = 0;

class point
{
public:
    double x;
    double y;
    double z;
    double w;

    point()
    {
        x = 0.0;
        y = 0.0;
        z = 0.0;
        w = 1.0;
    }

    point(double x1, double y1, double z1)
    {
        x = x1;
        y = y1;
        z = z1;
        w = 1.0;
    }

    point(double x1, double y1, double z1, double w1)
    {
        x = x1;
        y = y1;
        z = z1;
        w = w1;
    }

    point add_func(point my_point)
    {
        point temp;
        temp.x = x + my_point.x;
        temp.y = y + my_point.y;
        temp.z = z + my_point.z;
        return temp;
    }

    point sub_func(point my_point)
    {
        point temp;
        temp.x = x - my_point.x;
        temp.y = y - my_point.y;
        temp.z = z - my_point.z;
        return temp;
    }

    point mul_func(point my_point)
    {
        point temp;
        temp.x = x * my_point.x;
        temp.y = y * my_point.y;
        temp.z = z * my_point.z;
        return temp;
    }

    point div_func(point my_point)
    {
        point temp;
        temp.x = x / my_point.x;
        temp.y = y / my_point.y;
        temp.z = z / my_point.z;
        return temp;
    }

    point scalar_mul(double val)
    {
        point temp;
        temp.x = x * val;
        temp.y = y * val;
        temp.z = z * val;
        return temp;
    }

    double point_dot(point my_point)
    {
        double dot_val = 0.0;
        dot_val =  x*my_point.x + y*my_point.y + z*my_point.z;
        return dot_val;
    }

    point cross_product(point my_point)
    {
        point temp;
        temp.x = y*my_point.z - z*my_point.y;
        temp.y = z*my_point.x - x*my_point.z;
        temp.z = x*my_point.y - y*my_point.x;
        return temp;
    }

    void equal_func(point my_point)
    {
        x = my_point.x;
        y = my_point.y;
        z = my_point.z;
        w = my_point.w;
    }

    void scalePoint()
    {
        x = x*1.0/w;
        y = y*1.0/w;
        z = z*1.0/w;
        w = w*1.0/w;
    }

    void normalize()
    {
        double val;
        val = sqrt(x*x + y*y + z*z);
        x = x/val;
        y = y/val;
        z = z/val;
    }

    void printPoint(fstream& outfile)
    {
        cout.precision(7);
        outfile.precision(7);
        //cout<<"x is :"<<fixed<<x<<endl;
        //cout<<"y is :"<<fixed<<y<<endl;
        //cout<<"z is :"<<fixed<<z<<endl;
        outfile<<fixed<<x<<" ";
        outfile<<fixed<<y<<" ";
        outfile<<fixed<<z<<endl;
    }

    void printPoint()
    {
        cout.precision(7);
        cout<<"x is :"<<fixed<<x<<endl;
        cout<<"y is :"<<fixed<<y<<endl;
        cout<<"z is :"<<fixed<<z<<endl;
        cout<<"w is :"<<fixed<<w<<endl;
    }

    ~point()
    {
        x = 0.0;
        y = 0.0;
        z = 0.0;
        w = 1.0;
    }
};

class Matrix
{
public:
    double m[4][4];

    Matrix(int row, int col)
    {
        for(int i = 0; i < row; i++)
        {
            for(int j = 0; j < col; j++)
            {
                if(i == j)
                {
                    m[i][j] = 1;
                }
                else
                {
                    m[i][j] = 0;
                }
            }
        }
    }

    point Rodrigues(point x, point a, double angle)
    {
        point temp_1;
        point temp_2;
        point temp_3;
        point temp_4;

        temp_1.equal_func(x.scalar_mul(cos(angle*pi/180)));
        temp_2.equal_func(a.scalar_mul(a.point_dot(x)* (1 - cos(angle*pi/180))));
        temp_3.equal_func(a.cross_product(x).scalar_mul(sin(angle*pi/180)));

        temp_4.equal_func(temp_1.add_func(temp_2));
        temp_4.equal_func(temp_4.add_func(temp_3));

        return temp_4;
    }

    void genTranslationMatrix(double tx, double ty, double tz)
    {
        m[0][3] = tx;
        m[1][3] = ty;
        m[2][3] = tz;
    }

    void genScaleMatrix(double sx, double sy, double sz)
    {
        m[0][0] = sx;
        m[1][1] = sy;
        m[2][2] = sz;
    }

    void genRotateMatrix(double angle, double ax, double ay, double az)
    {
        point a(ax,ay,az);
        a.normalize();

        point i(1.0,0,0);
        point j(0,1.0,0);
        point k(0,0,1.0);

        point c1;
        point c2;
        point c3;

        c1.equal_func(Rodrigues(i,a,angle));
        c2.equal_func(Rodrigues(j,a,angle));
        c3.equal_func(Rodrigues(k,a,angle));

        m[0][0] = c1.x;
        m[1][0] = c1.y;
        m[2][0] = c1.z;

        m[0][1] = c2.x;
        m[1][1] = c2.y;
        m[2][1] = c2.z;

        m[0][2] = c3.x;
        m[1][2] = c3.y;
        m[2][2] = c3.z;
    }

    Matrix genViewMatrix(point eye, point look, point up)
    {
        point l;
        point r;
        point u;

        l.equal_func(look.sub_func(eye));
        l.normalize();
        r.equal_func(l.cross_product(up));
        r.normalize();
        u.equal_func(r.cross_product(l));

        Matrix T(4,4);
        Matrix R(4,4);
        Matrix V(4,4);

        T.m[0][3] = -eye.x;
        T.m[1][3] = -eye.y;
        T.m[2][3] = -eye.z;

        R.m[0][0] = r.x;
        R.m[0][1] = r.y;
        R.m[0][2] = r.z;
        R.m[1][0] = u.x;
        R.m[1][1] = u.y;
        R.m[1][2] = u.z;
        R.m[2][0] = -l.x;
        R.m[2][1] = -l.y;
        R.m[2][2] = -l.z;

        //T.printMatrix();
        //R.printMatrix();
        V.Mat_equal(mat_mul(R,T));
        return V;
    }

    void genProjectionMatrix(double fovY, double aspectratio, double near, double far)
    {
        double fovX,t,r;

        fovX = fovY*aspectratio;
        t = near*tan((fovY*1.0/2.0)*pi/180);
        r = near*tan((fovX*1.0/2.0)*pi/180);

        //point i(1.0,0,0);
        //point j(0,1.0,0);
        //point k(0,0,1.0);

        m[0][0] = near/r;
        m[1][1] = near/t;
        m[2][2] = -((far+near)/(far-near));
        m[2][3] = -((2*far*near)/(far-near));
        m[3][2] = -1;
        m[3][3] = 0;
    }

    point mul_point(point p)
    {
        double m1[4][1];
        double pnt[4][1];

        pnt[0][0] = p.x;
        pnt[1][0] = p.y;
        pnt[2][0] = p.z;
        pnt[3][0] = p.w;

        for (int i = 0; i < 4; i++)
        {
            m1[i][0] = 0.0;

            for (int j = 0; j < 4; j++)
            {
                m1[i][0] = m1[i][0] + (m[i][j] * pnt[j][0]);
            }
        }

        point temp(m1[0][0],m1[1][0],m1[2][0],m1[3][0]);
        return temp;
    }

    Matrix mat_mul(Matrix m1, Matrix m2)
    {
        Matrix temp(4,4);
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                temp.m[i][j] = 0.0;
                for (int k = 0; k < 4; k++)
                {
                    temp.m[i][j] = temp.m[i][j] + m1.m[i][k] * m2.m[k][j];
                }
            }
        }
        return temp;
    }

    void Mat_equal(Matrix temp)
    {
        for (int i = 0; i < 4 ; i++)
        {
            for (int j = 0 ; j < 4; j++)
            {
                m[i][j] = temp.m[i][j];
            }
        }
    }

    void printMatrix()
    {
        cout<<"matrix"<<endl<<endl;
        cout<<m[0][0]<<" "<<m[0][1]<<" "<<m[0][2]<<" "<<m[0][3]<<endl;
        cout<<m[1][0]<<" "<<m[1][1]<<" "<<m[1][2]<<" "<<m[1][3]<<endl;
        cout<<m[2][0]<<" "<<m[2][1]<<" "<<m[2][2]<<" "<<m[2][3]<<endl;
        cout<<m[3][0]<<" "<<m[3][1]<<" "<<m[3][2]<<" "<<m[3][3]<<endl;
    }

    ~Matrix()
    {
        for (int i = 0; i < 4; i++)
        {
            for (int j = 0; j < 4; j++)
            {
                m[i][j] = 0.0;
            }
        }
    }
};

class Color
{
public:
    int r;
    int g;
    int b;

    Color()
    {
        r = 0;
        g = 0;
        b = 0;
    }

    ~Color()
    {
        r = 0;
        g = 0;
        b = 0;
    }
};

class Triangle
{
public:
    point points[3];
    Color colors;

    Triangle()
    {
        for (int i = 0; i < 3; i++)
        {
            points[i].x = 0.0;
            points[i].y = 0.0;
            points[i].z = 0.0;
            points[i].w = 1.0;
        }
        colors.r = 0;
        colors.g = 0;
        colors.b = 0;
    }
};

void stage1(fstream& file, fstream& outfile_1)
{
    string command;
    stack<Matrix> s;
    stack<Matrix> pushmatrices;
    s.push(Matrix(4,4));

    while(true)
    {
        file>>command;

        if(command == "triangle")
        {
            point p1,p2,p3;

            file>>p1.x;
            file>>p1.y;
            file>>p1.z;

            file>>p2.x;
            file>>p2.y;
            file>>p2.z;

            file>>p3.x;
            file>>p3.y;
            file>>p3.z;

            p1.equal_func(s.top().mul_point(p1));
            p1.scalePoint();
            p2.equal_func(s.top().mul_point(p2));
            p2.scalePoint();
            p3.equal_func(s.top().mul_point(p3));
            p3.scalePoint();

            p1.printPoint(outfile_1);
            p2.printPoint(outfile_1);
            p3.printPoint(outfile_1);
            outfile_1<<endl;
            triangle++;
        }
        else if (command == "translate")
        {
            double tx,ty,tz;
            file>>tx;
            file>>ty;
            file>>tz;
            //cout<<tx<<endl<<ty<<endl<<tz<<endl;

            Matrix translationMat(4,4);
            translationMat.genTranslationMatrix(tx,ty,tz);

            Matrix temp(4,4);
            temp.Mat_equal(temp.mat_mul(s.top(),translationMat));
            s.pop();
            s.push(temp);
        }
        else if (command == "scale")
        {
            double sx,sy,sz;
            file>>sx;
            file>>sy;
            file>>sz;
            //cout<<sx<<endl<<sy<<endl<<sz<<endl;

            Matrix scaleMat(4,4);
            scaleMat.genScaleMatrix(sx,sy,sz);

            Matrix temp(4,4);
            temp.Mat_equal(temp.mat_mul(s.top(),scaleMat));
            s.pop();
            s.push(temp);
        }
        else if (command == "rotate")
        {
            double angle;
            double ax,ay,az;
            file>>angle;
            file>>ax;
            file>>ay;
            file>>az;
            //cout<<angle<<endl;
            //cout<<ax<<endl<<ay<<endl<<az<<endl;

            Matrix rotateMat(4,4);
            rotateMat.genRotateMatrix(angle,ax,ay,az);
            //rotateMat.printMatrix();

            Matrix temp(4,4);
            //s.top().printMatrix();
            temp.Mat_equal(temp.mat_mul(s.top(),rotateMat));
            //temp.printMatrix();
            s.pop();
            s.push(temp);
        }
        else if (command == "push")
        {
            //cout<<"push"<<endl;
            /*Matrix temp(4,4);
            temp.Mat_equal(temp.mat_mul(temp,s.top()));
            s.push(temp);*/
            pushmatrices.push(s.top());
            push_count++;
        }
        else if (command == "pop")
        {
            //cout<<"pop"<<endl;
            if(push_count == 0)
            {
                cout<<"Pop matrix called before push function"<<endl;
                exit(1);
            }
            s.pop();
            s.push(pushmatrices.top());
            pushmatrices.pop();
            push_count--;
        }
        else if (command == "end")
        {
            //cout<<"end"<<endl;
            break;
        }
    }
    file.close();
    outfile_1.close();
}

void stage2(point eye, point look, point up, fstream& file, fstream& outfile_2)
{
    Matrix viewMatrix(4,4);
    viewMatrix.Mat_equal(viewMatrix.genViewMatrix(eye,look,up));

    file.open("stage1.txt",ios::in);
    if(!file.is_open()){
		cout<<"Cannot open file"<<endl;
        exit(1);
	}

    //cout<<"stage_2"<<endl;
    //viewMatrix.printMatrix();
    for (int i = 0; i < triangle; i++)
    {
        point p1,p2,p3;

        file>>p1.x;
        file>>p1.y;
        file>>p1.z;

        file>>p2.x;
        file>>p2.y;
        file>>p2.z;

        file>>p3.x;
        file>>p3.y;
        file>>p3.z;

        p1.equal_func(viewMatrix.mul_point(p1));
        p1.scalePoint();
        p2.equal_func(viewMatrix.mul_point(p2));
        p2.scalePoint();
        p3.equal_func(viewMatrix.mul_point(p3));
        p3.scalePoint();

        p1.printPoint(outfile_2);
        p2.printPoint(outfile_2);
        p3.printPoint(outfile_2);
        outfile_2<<endl;
    }
    file.close();
    outfile_2.close();
}

void stage3(double fovY, double aspectRatio, double near, double far, fstream& file, fstream& outfile_3)
{
    Matrix projectMatrix(4,4);
    projectMatrix.genProjectionMatrix(fovY,aspectRatio,near,far);
    //projectMatrix.printMatrix();

    file.open("stage2.txt",ios::in);
    if(!file.is_open()){
		cout<<"Cannot open file"<<endl;
        exit(1);
	}

    //cout<<"stage_3"<<endl;
    //projectMatrix.printMatrix();
    for (int i = 0; i < triangle; i++)
    {
        point p1,p2,p3;

        file>>p1.x;
        file>>p1.y;
        file>>p1.z;

        file>>p2.x;
        file>>p2.y;
        file>>p2.z;

        file>>p3.x;
        file>>p3.y;
        file>>p3.z;

        p1.equal_func(projectMatrix.mul_point(p1));
        //p1.printPoint();
        //cout<<endl;
        p1.scalePoint();
        p2.equal_func(projectMatrix.mul_point(p2));
        p2.scalePoint();
        p3.equal_func(projectMatrix.mul_point(p3));
        p3.scalePoint();

        p1.printPoint(outfile_3);
        p2.printPoint(outfile_3);
        p3.printPoint(outfile_3);
        outfile_3<<endl;
    }
    file.close();
    outfile_3.close();
}

double max_val(double x, double y, double z)
{
    double temp;
    temp = max(x,y);
    temp = max(temp,z);
    return temp;
}

double min_val(double x, double y, double z)
{
    double temp;
    temp = min(x,y);
    temp = min(temp,z);
    return temp;
}

int convertToint(double val)
{
    return (int)round(val);
}

int get_topScaline(double TopY, double max_Y, double dy)
{
    int top = 0;
    if  (max_Y < TopY)
    {
        //cout<<"top : "<<TopY<<endl<<"max "<<max_Y<<endl<<"dy: "<<dy<<endl;
        double pixel_num = (TopY - max_Y)/dy;
        top = convertToint(pixel_num);
    }
    return top;
}

int get_bottomScaline(double screenHeight, double BottomY, double min_Y, double dy)
{
    int bottom = screenHeight - 1;
    if (min_Y > BottomY)
    {
        double pixel_num = (min_Y - BottomY)*1.0/dy;
        bottom = screenHeight - (1 + convertToint(pixel_num));
    }
    return bottom;
}

int getLeftIntersectingColumn(double LeftX, double min_X, double dx)
{
    int left_val = 0;
    if (min_X > LeftX)
    {
        double pixel_num = (min_X - LeftX)*1.0/dx;
        left_val = convertToint(pixel_num);
    }
    return left_val;
}

int getRightIntersectingColumn(double screenWidth, double RightX, double max_X, double dx)
{
    int right_val = screenWidth - 1;
    if (max_X < RightX)
    {
        double pixel_num = (RightX - max_X)*1.0/dx;
        right_val = screenWidth - (1 + convertToint(pixel_num));
    }
    return right_val;
}

int getMinI(point intersections[3], int minI)
{
    double min_X;

    int i = 0;
    while(i < 3)
    {
        if (minI == -1)
        {
            if (intersections[i].x != inf)
            {
                minI = i;
                min_X = intersections[i].x;
            }
        }
        else
        {
            if (intersections[i].x != inf)
            {
                if (intersections[i].x < min_X)
                {
                    minI = i;
                    min_X = intersections[i].x;
                }
            }
        }
        i++;
    }
    return minI;
}

int getMaxI(point intersections[3], int maxI)
{
    double max_X;

    int i = 0;
    while(i < 3)
    {
        if (maxI == -1)
        {
            if (intersections[i].x != inf)
            {
                maxI = i;
                max_X = intersections[i].x;
            }
        }
        else
        {
            if (intersections[i].x != inf)
            {
                if(intersections[i].x > max_X)
                {
                    maxI = i;
                    max_X = intersections[i].x;
                }
            }
        }
        i++;
    }
    return maxI;
}

void stage4(fstream& file)
{
    fstream outfile_4;
    int screenWidth,screenHeight;
    double leftLimitX,rightLimitX,topLimitY,bottomLimitY,frontLimitZ,rearLimitZ;

    file.open("config.txt",ios::in);
    if(!file.is_open()){
		cout<<"Cannot open file"<<endl;
        exit(1);
	}

	file>>screenWidth;
	file>>screenHeight;
	file>>leftLimitX;
	file>>bottomLimitY;
	file>>frontLimitZ;
	file>>rearLimitZ;

	//cout<<screenWidth<<endl<<screenHeight<<endl<<leftLimitX<<endl<<bottomLimitY<<endl<<frontLimitZ<<endl<<rearLimitZ<<endl;

	rightLimitX = -leftLimitX;
	topLimitY = -bottomLimitY;

	file.close();

	Triangle my_triangles[triangle];

	file.open("stage3.txt",ios::in);
    if(!file.is_open()){
		cout<<"Cannot open file"<<endl;
        exit(1);
	}

    srand(time(NULL));
	for (int i = 0; i < triangle; i++)
    {
        file>>my_triangles[i].points[0].x;
        file>>my_triangles[i].points[0].y;
        file>>my_triangles[i].points[0].z;
        file>>my_triangles[i].points[1].x;
        file>>my_triangles[i].points[1].y;
        file>>my_triangles[i].points[1].z;
        file>>my_triangles[i].points[2].x;
        file>>my_triangles[i].points[2].y;
        file>>my_triangles[i].points[2].z;

        //my_triangles[i].points[0].printPoint();
        //my_triangles[i].points[1].printPoint();
        //my_triangles[i].points[2].printPoint();

        my_triangles[i].colors.r = rand()%256;
        my_triangles[i].colors.g = rand()%256;
        my_triangles[i].colors.b = rand()%256;
    }
    file.close();

    double dx = (rightLimitX - leftLimitX)/screenWidth;
    double dy = (topLimitY - bottomLimitY)/screenHeight;
    double TopY = topLimitY - (dy*1.0/2.0);
    double LeftX = leftLimitX + (dx*1.0/2.0);
    double BottomY = bottomLimitY + (dy*1.0/2.0);
    double RightX = rightLimitX - (dx*1.0/2.0);

    double** z_Buffer;
    z_Buffer = new double*[screenHeight];
    for (int i = 0; i < screenHeight; i++)
    {
        z_Buffer[i] = new double[screenWidth];
    }

    for (int i = 0; i < screenHeight; i++)
    {
        for (int j = 0; j < screenWidth; j++)
        {
            z_Buffer[i][j] = rearLimitZ;
        }
    }

    Color** frames;
    frames = new Color*[screenHeight];
    for (int i = 0; i < screenHeight; i++)
    {
        frames[i] = new Color[screenWidth];
    }

    for (int i = 0; i < screenHeight; i++)
    {
        for (int j = 0; j < screenWidth; j++)
        {
            frames[i][j].r = 0;
            frames[i][j].g = 0;
            frames[i][j].b = 0;
        }
    }

    int top_Scanline = 0;
    int bottom_Scanline = 0;


    for (int i = 0; i < triangle; i++)
    {
        double max_Y,min_Y;
        max_Y = max_val(my_triangles[i].points[0].y, my_triangles[i].points[1].y, my_triangles[i].points[2].y);
        min_Y = min_val(my_triangles[i].points[0].y, my_triangles[i].points[1].y, my_triangles[i].points[2].y);

        //cout<<"max : "<<max_Y<<endl<<"min : "<<min_Y<<endl;
        top_Scanline = get_topScaline(TopY,max_Y,dy);
        bottom_Scanline = get_bottomScaline(screenHeight,BottomY,min_Y,dy);

        int left_intersecting_column,right_intersecting_column;
        left_intersecting_column = 0;
        right_intersecting_column = 0;

        double y_val;
        //cout<<"top : "<<top_Scanline<<endl<<"bottom : "<<bottom_Scanline<<endl;
        for (int i1 = top_Scanline; i1 <= bottom_Scanline; i1++)
        {
            y_val = TopY - i1*dy;

            point intersections[3];
            for (int inpi = 0; inpi <= 2; inpi++)
            {
                intersections[inpi] = point(inf, y_val, inpi, inpi+1);
            }
            intersections[2].w = 0;


            for (int j = 0; j < 3; j++)
            {
                point p1;
                p1.equal_func(my_triangles[i].points[(int)intersections[j].z]);
                point p2;
                p2.equal_func(my_triangles[i].points[(int)intersections[j].w]);

                if (p1.y != p2.y)
                {
                    intersections[j].x = (p1.x + (y_val - p1.y)*(p1.x - p2.x)/(p1.y - p2.y));
                }
            }

            for (int j = 0; j < 3; j++)
            {
                point p1;
                p1.equal_func(my_triangles[i].points[(int)intersections[j].z]);
                point p2;
                p2.equal_func(my_triangles[i].points[(int)intersections[j].w]);

                if (intersections[j].x != inf)
                {
                    bool cond_1 = intersections[j].x > max(p1.x,p2.x);
                    bool cond_2 = intersections[j].x < min(p1.x,p2.x);
                    bool cond_3 = intersections[j].y > max(p1.y,p2.y);
                    bool cond_4 = intersections[j].y < min(p1.y,p2.y);
                    if (cond_1 || cond_2 || cond_3 || cond_4)
                    {
                        intersections[j].x = inf;
                    }
                }
            }


            int maxI;
            int minI;

            minI = getMinI(intersections,-1);
            maxI = getMaxI(intersections,-1);


            //cout<<"max_x : "<<max_X<<endl<<"inter : "<<intersections[maxI].x<<endl;
            //cout<<"min_x : "<<min_X<<endl<<"inter : "<<intersections[minI].x<<endl;
            left_intersecting_column = getLeftIntersectingColumn(LeftX,intersections[minI].x,dx);
            right_intersecting_column = getRightIntersectingColumn(screenWidth,RightX,intersections[maxI].x,dx);

            point p1;
            p1.equal_func(my_triangles[i].points[(int)intersections[minI].z]);
            point p2;
            p2.equal_func(my_triangles[i].points[(int)intersections[minI].w]);

            double za;
            za = p1.z + (intersections[minI].y - p1.y)*(p2.z - p1.z)/(p2.y - p1.y);

            p1.equal_func(my_triangles[i].points[(int)intersections[maxI].z]);
            p2.equal_func(my_triangles[i].points[(int)intersections[maxI].w]);

            double zb;
            zb = p1.z + (intersections[maxI].y - p1.y)*(p2.z - p1.z)/(p2.y - p1.y);

            double zp;
            double xb = intersections[maxI].x;
            double xa = intersections[minI].x;
            double cons_val;
            cons_val = dx*(zb - za)/(xb - xa);

            //cout<<"left : "<<left_intersecting_column<<endl<<"right : "<<right_intersecting_column<<endl;
            //leftfile<<"left : "<<left_intersecting_column<<endl<<"right : "<<right_intersecting_column<<endl;
            for (int j1 = left_intersecting_column; j1 <= right_intersecting_column; j1++)
            {

                if (j1 == left_intersecting_column)
                {
                    zp = za + ((LeftX + left_intersecting_column*dx) - xa)*(zb - za)/(xb - xa);
                }
                else
                {
                    zp = zp + cons_val;
                }

                if ( zp > frontLimitZ && zp < z_Buffer[i1][j1])
                {
                    z_Buffer[i1][j1] = zp;

                    frames[i1][j1].r = my_triangles[i].colors.r;
                    frames[i1][j1].g = my_triangles[i].colors.g;
                    frames[i1][j1].b = my_triangles[i].colors.b;
                }
            }
        }
    }

    bitmap_image bitmapImage(screenWidth, screenHeight);

    for (int i = 0; i < screenHeight; i++)
    {
        for (int j = 0; j < screenWidth; j++)
        {
            bitmapImage.set_pixel(j, i, frames[i][j].r, frames[i][j].g, frames[i][j].b);
        }
    }

    bitmapImage.save_image("out.bmp");

    outfile_4.open("z_buffer.txt",ios::out);
    if(!outfile_4.is_open()){
		cout<<"Cannot open file"<<endl;
        exit(1);
	}

    outfile_4.precision(6);
    for (int i = 0; i < screenHeight; i++)
    {
        for (int j = 0; j < screenWidth; j++)
        {
            if (z_Buffer[i][j] < rearLimitZ)
            {
                outfile_4<<fixed<<z_Buffer[i][j];
                outfile_4<<'\t';
            }
        }
        outfile_4<<endl;
    }
    outfile_4.close();

    for (int i = 0; i < screenHeight; i++)
    {
        delete[] z_Buffer[i];
    }
    delete[] z_Buffer;

    for (int i = 0; i < screenHeight; i++)
    {
        delete[] frames[i];
    }
    delete[] frames;

}

int main(int argc, char *argv[])
{
    //cout<<INF<<endl;
    fstream file;
    fstream outfile_1;
    fstream outfile_2;
    fstream outfile_3;
    file.open("scene.txt",ios::in);
	if(!file.is_open()){
		cout<<"Cannot open file"<<endl;
        exit(1);
	}

	outfile_1.open("stage1.txt",ios::out);
	if(!outfile_1.is_open()){
		cout<<"Cannot open file"<<endl;
        exit(1);
	}

	outfile_2.open("stage2.txt",ios::out);
	if(!outfile_2.is_open()){
		cout<<"Cannot open file"<<endl;
        exit(1);
	}

	outfile_3.open("stage3.txt",ios::out);
	if(!outfile_3.is_open()){
		cout<<"Cannot open file"<<endl;
        exit(1);
	}

    double eyeX,eyeY,eyeZ;
    double lookX,lookY,lookZ;
    double upX,upY,upZ;
    double fovY,aspectRatio,near,far;

    file>>eyeX;
    file>>eyeY;
    file>>eyeZ;

    point eye(eyeX,eyeY,eyeZ);
    //eye.printPoint();

    file>>lookX;
    file>>lookY;
    file>>lookZ;

    point look(lookX,lookY,lookZ);
    //look.printPoint();

    file>>upX;
    file>>upY;
    file>>upZ;

    point up(upX,upY,upZ);
    //up.printPoint();

    file>>fovY;
    file>>aspectRatio;
    file>>near;
    file>>far;

    //cout<<fovY<<endl<<aspectRatio<<endl<<near<<endl<<far<<endl;

    stage1(file,outfile_1);
    stage2(eye,look,up,file,outfile_2);
    stage3(fovY,aspectRatio,near,far,file,outfile_3);
    stage4(file);

	return 0;
}

