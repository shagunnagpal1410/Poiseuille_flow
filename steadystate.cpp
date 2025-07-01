#include<iostream>
#include<vector>
#include<map>
#include<cmath>
#include<Eigen/Sparse>
#include<Eigen/Dense>
#include<Eigen/SparseLU>
#include<fstream>
#include<string>
using namespace std;
using namespace Eigen;
typedef Triplet<double> Tri;
struct point{
    double x,y;
    int voxel;
    point(double i, double j) {
        x=i, y=j;
        voxel=0;
    }
    bool operator<(const point& other) const {
        if (x != other.x) return x < other.x;
        if (y != other.y) return y < other.y;
        return voxel < other.voxel;
    }
};
int find_voxel(double x, double y, double x_min, double x_max, int voxels_inrow, double radius) {
    return int((x-x_min)/radius) + int((y)/radius)*voxels_inrow;
}
int find_rownumber(int voxel_number, int voxels_inrow) {
    return voxel_number%voxels_inrow;
}
int find_columnnumber(int voxel_number, int voxels_inrow) {
    return int(voxel_number/voxels_inrow);
}
double gaussian_weight_function(point Ni, point p0, double radius) {
    double distance=(pow(Ni.x-p0.x,2)+pow(Ni.y-p0.y,2))/(radius*radius);
    if (distance<=1) {
        return exp(-6.25*distance);
    }
    else {
        return 0.0;
    }
}
int main() {
    //defining domain
    double L=2.0;
    double B=1.0;
    double dx=0.1;
    double dy=0.1;
    vector<point> previous_domain;
    double x_min=0.0;
    double x_max=2.0;
    double y_min=0.0;
    double y_max=1.0;
    double radius=0.3;
    int voxels_inrow=int((x_max-x_min)/radius)+1;
    int voxels_incolumn=int(1/radius)+1;
    //initializing some parameters------------------Ending---------------------------------------------------------
    //adding points in the domain-------------------Starting-------------------------------------------------------
    int n=ceil((x_max-x_min)/dx);
    int m=ceil((y_max-y_min)/dy);
    for (int i=0; i<=n; i++) {
        for (int j=0; j<=m; j++) {
            point p1=point(i*dx,j*dy);
            p1.voxel=find_voxel(p1.x,p1.y,x_min,x_max,voxels_inrow,radius);
            previous_domain.push_back(p1);
        }
    }
    int max_voxels=voxels_incolumn*voxels_inrow;
    vector<vector<point>> points_insidevoxel(max_voxels);
    for (int i=0; i<previous_domain.size(); i++) {
        point p0=previous_domain[i];
        points_insidevoxel[p0.voxel].push_back(p0);
    }
    //adding points in the domain-------------------Ending---------------------------------------------------------
    //finding neighbour voxels of each voxel--------Starting-------------------------------------------------------
    vector<vector<int>> neighbour_voxels(max_voxels);
    for(int i=0; i<max_voxels; i++) {
        int find_x=find_rownumber(i, voxels_inrow);
        int find_y=find_columnnumber(i, voxels_inrow);
        for (int diffx=-1; diffx<=1; diffx++) {
            for (int diffy=-1; diffy<=1; diffy++) {
                if (find_x+diffx>=0 && find_x+diffx<=voxels_inrow-1 && diffy+find_y>=0 && find_y+diffy<=voxels_incolumn-1) {
                            neighbour_voxels[i].push_back((find_y+diffy)*voxels_inrow+(find_x+diffx));
                }
            }
        }
    }
    //finding neighbour voxels of each voxel--------Ending--------------------------------------------------------
    //we will be finding neighbours of each point respectively
    map<point,vector<point>> neighbours_ofpoint;
    for (int p=0; p<previous_domain.size(); p++) {
        point p0=previous_domain[p];
        int voxel_number=p0.voxel;
        for (int i=0; i<neighbour_voxels[voxel_number].size(); i++) {
            int neighbour_voxel=neighbour_voxels[voxel_number][i];
            for (int j=0; j<points_insidevoxel[neighbour_voxel].size(); j++) {
                point Ni=points_insidevoxel[neighbour_voxel][j];
                double distance=pow(Ni.x-p0.x,2)+pow(Ni.y-p0.y,2);
                if (distance>0 && distance<pow(radius,2)) {
                    neighbours_ofpoint[p0].push_back(Ni);
                }
            }
        }
    }
    map<point,double> velocity;
        int total_points=previous_domain.size();
        SparseMatrix<double> letsSolve(total_points,total_points);
        VectorXd rhs(total_points);
        vector<Tri> coefficients;
        map<point,int> identity;
        for (int p=0; p<previous_domain.size(); p++) {
            identity[previous_domain[p]]=p;
        }
        for (int p=0; p<previous_domain.size(); p++) {
            point p0=previous_domain[p];
            vector<point> neighbours=neighbours_ofpoint[p0];
            int total_neighbours=neighbours.size();
            if (p0.y==y_min || p0.y==y_max) {
                MatrixXd M(total_neighbours+2,6);
                MatrixXd W=MatrixXd :: Zero(total_neighbours+2,total_neighbours+2);
                for (int i=0; i<total_neighbours; i++) {
                    point Ni=neighbours[i];
                    double Dx=Ni.x-p0.x, Dy=Ni.y-p0.y;
                    M(i,0)=1, M(i,1)=Dx, M(i,2)=Dy, M(i,3)=(Dx*Dx)/2.0, M(i,4)=(Dy*Dy)/2.0, M(i,5)=(Dx*Dy);
                    W(i,i)=gaussian_weight_function(Ni,p0,radius);
                }
                M(total_neighbours,0)=0,M(total_neighbours,1)=0,M(total_neighbours,2)=0,M(total_neighbours,3)=1
                ,M(total_neighbours,4)=1,M(total_neighbours,5)=0; 
                W(total_neighbours,total_neighbours)=1;
                M(total_neighbours+1,0)=1,M(total_neighbours+1,1)=0,M(total_neighbours+1,2)=0,M(total_neighbours+1,3)=0
                ,M(total_neighbours+1,4)=0,M(total_neighbours+1,5)=0; 
                W(total_neighbours+1,total_neighbours+1)=1;
                MatrixXd MTWM=M.transpose()*W*M;
                MatrixXd MTW=M.transpose()*W;  
                MatrixXd A=MTWM.ldlt().solve(MTW);
                double constant=0.0;
                for (int i=0; i<total_neighbours; i++) {
                    point Ni=neighbours[i];
                    coefficients.push_back(Tri(identity[p0],identity[Ni], A(0,i)));
                }
                constant=constant-(A(0,total_neighbours)*1.0);
                constant=constant-(A(0,total_neighbours+1)*0.0);
                coefficients.push_back(Tri(identity[p0],identity[p0],-1));
                rhs(identity[p0])=constant;
            }
            else {
                MatrixXd M(total_neighbours+1,6);
                MatrixXd W=MatrixXd :: Zero(total_neighbours+1,total_neighbours+1);
                for (int i=0; i<total_neighbours; i++) {
                    point Ni=neighbours[i];
                    double Dx=Ni.x-p0.x, Dy=Ni.y-p0.y;
                    M(i,0)=1, M(i,1)=Dx, M(i,2)=Dy, M(i,3)=(Dx*Dx)/2.0, M(i,4)=(Dy*Dy)/2.0, M(i,5)=(Dx*Dy);
                    W(i,i)=gaussian_weight_function(Ni,p0,radius);
                }
                M(total_neighbours,0)=0,M(total_neighbours,1)=0,M(total_neighbours,2)=0,M(total_neighbours,3)=1
                ,M(total_neighbours,4)=1,M(total_neighbours,5)=0; 
                W(total_neighbours,total_neighbours)=1;
                MatrixXd MTWM=M.transpose()*W*M;
                MatrixXd MTW=M.transpose()*W;  
                MatrixXd A=MTWM.ldlt().solve(MTW);
                double constant=0.0;
                for (int i=0; i<total_neighbours; i++) {
                    point Ni=neighbours[i];
                    coefficients.push_back(Tri(identity[p0],identity[Ni], A(0,i)));
                }
                constant=constant-(A(0,total_neighbours)*1.0);
                coefficients.push_back(Tri(identity[p0],identity[p0],-1));
                rhs(identity[p0])=constant;
            }
        }
        SparseLU<SparseMatrix<double>, COLAMDOrdering<int>> solver;
        letsSolve.setFromTriplets(coefficients.begin(), coefficients.end());
        solver.compute(letsSolve);
        VectorXd answer = solver.solve(rhs);
        for (int p=0; p<previous_domain.size(); p++) {
            point p0=previous_domain[p];
            velocity[p0]=answer(identity[p0]);
        }
        ofstream fout("Velocity.csv");
        fout<<"X"<<","<<"Y"<<","<<"Velocity"<<"\n";
        for (int p=0; p<previous_domain.size(); p++) {
            point p0=previous_domain[p];
            fout<<p0.x<<","<<p0.y<<","<<-1*velocity[p0]<<endl;
        }
        fout.close();
    return 0;
}
