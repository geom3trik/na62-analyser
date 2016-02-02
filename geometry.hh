#ifndef GEOMETRY_HH_INCLUDED
#define GEOMETRY_HH_INCLUDED

struct Intersect
{
    bool flag;
    TVector3 first; //First intersection point
    TVector3 second; //Second intersection point
};

double sq(double value)
{
    return value*value;
}

class Ray3
{
    public:
        Ray3(TVector3 origin_, TVector3 direction_) : ori(origin_), dir(direction_){}
        TVector3 ori;
        TVector3 dir;

};

class Box3
{
    public:
        Box3(TVector3 lo_, TVector3 hi_) : lo(lo_), hi(hi_){}
        TVector3 lo;
        TVector3 hi;
};

class Plane
{
    public:
        Plane(TVector3 origin_, TVector3 normal_) : origin(origin_), normal(normal_){}
        TVector3 origin;
        TVector3 normal;
};

class Quadratic
{
    public:
        Quadratic(double a_=0, double b_=0, double c_=0) : a(a_), b(b_), c(c_), complex(false), x1(0), x2(0)
        {
            double det = pow(b,2) - 4*a*c;
            if(det >= 0)
            {
                x1 = (2*b+sqrt(det))/(4*a);
                x1 = (2*b-sqrt(det))/(4*a);
            }
            else
            {
                complex = true;
            }
        }
        bool complex;
        double a, b, c;
        double x1;
        double x2;
};

class Cylinder
{
    public:
        Cylinder(TVector3 start_, TVector3 end_ , double radius_) : start(start_), end(end_), radius(radius_)
        {

        }
        Intersect intersect(Ray3& ray)
        {
            Intersect i;
            std::vector<double> results;
            TVector3 v = end - start;

            double a = (1/v.Mag2())*sq(ray.dir.Dot(v))-ray.dir.Mag2();
            double b = 2*(1/v.Mag2())*ray.ori.Dot(v)*ray.dir.Dot(v)-2*ray.ori.Dot(ray.dir);
            double c = sq(radius)+(1/v.Mag2())*sq(ray.ori.Dot(v));
            Quadratic t(a,b,c);
            //Only interested in positive real values
            if(!t.complex)
            {
                if(t.x1 >= 0)
                {
                    i.flag = true;
                    //i.first = ray.ori + t.x1*ray.dir;
                    results.push_back(t.x1);
                }
                if(t.x2 >= 0)
                {
                    i.flag = true;
                    //i.second = ray.ori + x2*ray.dir;
                    results.push_back(t.x2);
                }
            }
            double cap1 = (-v.Dot(ray.ori))/(v.Dot(ray.dir));
            TVector3 p = ray.ori + cap1*ray.dir;
            if(p.Dot(p) <= sq(radius))
            {
                results.push_back(cap1);
            }
            double cap2 = (v.Mag2()-v.Dot(ray.ori))/(v.Dot(ray.dir));
            p = ray.ori + cap2*ray.dir;
            TVector3 m = p - v;
            if(m.Mag2() <= sq(radius))
            {
                results.push_back(cap2);
            }
            if(results.size() == 1 || results.size() > 2 || results.size() < 0)
            {
                std::cout << "something went wrong!" << std::endl;
            }
            else if(results.size() == 2)
            {
                i.first = ray.ori + results.front()*ray.dir;
                i.second = ray.ori + results.back()*ray.dir;
                i.flag = true;
            }
            else
            {
                i.flag = false;
            }


        }
        TVector3 start;
        TVector3 end;
        double radius;
};

/*
class RayIntersect
{
    public:
        RayIntersect(Ray& ray)

        RayIntersect(Ray& ray, Box3& box)
        {
            double tmin = -INFINITY, tmax = INFINITY;

            if(ray.direction[0] != 0)
            {
                double tx1 = (box.lo[0] - ray.origin[0])/ray.direction[0];
                double tx2 = (box.lo[0] - ray.origin[0])/ray.direction[0];

                tmin = max(tmin,min(tx1,tx2));
                tmax = min(tmax,max(tx1,tx2));
            }
            if(ray.direction[0] != 0)
            {
                double ty1 = (box.lo[1] - ray.origin[1])/ray.direction[1];
                double ty2 = (box.lo[1] - ray.origin[1])/ray.direction[1];

                tmin = max(tmin,min(ty1,ty2));
                tmax = min(tmax,max(ty1,ty2));
            }
            if(ray.direction[0] != 0)
            {
                double tz1 = (box.lo[2] - ray.origin[2])/ray.direction[2];
                double tz2 = (box.lo[2] - ray.origin[2])/ray.direction[2];

                tmin = max(tmin,min(tx1,tx2));
                tmax = min(tmax,max(tx1,tx2));
            }

            intersected = (tmax >= tmin);
            length = sqrt(pow(tx2-tx1,2)+pow(ty2-ty1,2)+pow(tz2-tz1,2));
        }
        bool intersect() const {return intersected;}
        double length() const {return length;}
    private:
        bool intersect_;
        double length_;
};
*/



#endif // GEOMETRY_HH_INCLUDED
