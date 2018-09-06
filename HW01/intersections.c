#include <math.h>
#include "intersections.h"
#include <stdio.h>

//if there is no intersection, this function should return 0
//Otherwise, populate 'intersection' with the point of intersection and return 1
int ray_sphere_intersection(ray_t observer, sphere_t obj, vector_t *intersection) {
  
  //linear term inside the norm
  vector_t h_term = observer.dir;
  //constant term inside the norm
  vector_t const_term = difference(observer.start,obj.center);
  
  //quadratic term in the expansion
  double a = dot_product(h_term,h_term);
  //linear term in the expansion
  double b = 2*dot_product(h_term,const_term);
  //constant term in the expansion
  double c = dot_product(const_term,const_term) - obj.radius*obj.radius;

  double discriminant = b*b - 4*a*c;

  //make sure the equation has real roots
  if (discriminant < 0.) return 0;
  //find the smallest positive h
  double h;
  if ((-1*b - sqrt(discriminant))/(2*a) > 0) h = (-1*b - sqrt(discriminant))/(2*a); //sphere in front of us
  else if ((-1*b + sqrt(discriminant))/(2*a) > 0) h = (-1*b + sqrt(discriminant))/(2*a); //we are inside the sphere
  else return 0; //no positive roots, so the sphere is behind us

  vector_t solution = scaled_sum(1.,observer.start,h,observer.dir);

  copy_vector(solution,intersection);
  return 1;  
}

//if there is no intersection, this function should return 0
//Otherwise, populate 'intersection' with the point of intersection and return 1
int ray_disk_intersection(ray_t observer, disk_t obj, vector_t *intersection) {
  //Question 3: Modify this function to compute an intersection

  //d; direction term
  vector_t d = observer.dir;
  //s; start term
  vector_t s = observer.start;

  float t;

//if n dot d is 0, vector is tangent to disk - doesn't count
  if(dot_product(n,d)==0){
    return 0;
  }else{

    //calculate t using formula derived in latex doc
    t = (dot_product(n,difference(c-s)))/(dot_product(n,d));

    //if t is negative, vector points away from disk
    if(t<0){
      return 0;
    } else{

      //find p and check distance from center
      vector_t p = sum(s+scalar_product(t,d));
      double dist = distance(obj.center,p);

      if(dist<=obj.radius){
	copy_vector(p,intersection);
	return 1;
      }
      
    }

  }
  
  return 0;
}

//if there is no intersection, this function should return 0
//Otherwise, populate 'intersection' with the point of intersection and return 1
int ray_cylinder_intersection(ray_t observer, cylinder_t obj, vector_t *intersection) {
  //Question 5: Modify this function to compute an intersection

  //d; direction term
  vector_t d = observer.dir;
  //s; start term
  vector_t s = observer.start;

  //a; axis
  vector_t a = obj.axis;
  //c; center
  vector_t c = obj.center;
  //r; radius
  double r = obj.radius;
  //h; height
  double h = obj.height;

  //take discriminant
  //s_tilda is a substitute variable - holds value
  //s-c-a(a*s)+a(a*c)
  vector_t s_tilda = difference(s,difference(c,sum(scalar_product(dot_product(a,s),a),scalar_product(dot_product(a,c),a))));
  float dDotD = dot_product(d,d);
  float stDotSt = dot_product(s_tilda,s_tilda);
  float discriminant = 2*dot_product(d,dDotD)-(4*dot_product(dDotD,stDotSt));

  //check if disc is <0
  if (discriminant<0){
    return 0;
  }

  //find p
  vector_t p = sum(s,scalar_product(t,d));
  
  //is point of intersection p between base point c and c+h

  //find point on c base closest to p
  //drop vector starting at p parallel to axis
  ray_t tempRay = new_ray(p,-a);
  //find intersection point with base disk
  disk_t tempDisk = new_disk(r,c,a);
  //p_tilda will hold point on base closest to p
  vector_t p_tilda;

  ray_disk_intersection(tempRay,tempDisk,p_tilda);

  //now is p between p_tilda and h
  float dist = distance(p_tilda,p);
  if(dist<=h){
    copy_vector(p,intersection);
    return 1;
  }
  
  return 0;
}

int ray_cone_intersection(ray_t observer, cone_t obj, vector_t *intersection) {
  //Question 7: Modify this function to compute an intersection

    //d; direction term
  vector_t d = observer.dir;
  //s; start term
  vector_t s = observer.start;

  //a; axis
  vector_t a = obj.axis;
  //c; center
  vector_t c = obj.center;
  //r; radius
  double r = obj.radius;
  //h; height
  double h = obj.height;

  //take quadratic
  float div = r/h;
  div = div*div;
  float a_term = d.x*d.x + d.y*d.y - div*d.z*d.z;
  float b_term = 2*(d.x*s.x + d.y*s.y - div*d.z*(s.z-h));
  float c_term = s.x*s.x + s.y*s.y - div*(s.z-h)*(s.z-h);
  float discriminant = b_term*b_term - 4*a_term*c_term;

  //check if disc is <0
  if (discriminant<0){
    return 0;//no intersection
  }

  //find p
  vector_t p = sum(s,scalar_product(t,d));
  
  //is point of intersection p between dummy point b and b-h

  //find point on c base closest to p
  //drop vector starting at p parallel to axis
  ray_t tempRay = new_ray(p,-a);
  //create dummy disk normal to axis a with center at vertex v
  disk_t tempDisk = new_disk(r,v,a);
  //p_tilda will hold point on base closest to p
  vector_t p_tilda;

  ray_disk_intersection(tempRay,tempDisk,p_tilda);

  //now is p between p_tilda and h
  float dist = distance(p_tilda,p);
  if(dist<=h){
    copy_vector(p,intersection);
    return 1;
  }
  
  return 0;
}
  
  return 0;
}
