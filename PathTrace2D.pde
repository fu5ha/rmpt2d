int MAX_RAY_DEPTH = 5; // How many bounces each original light beam can have
float MAX_RAY_DIST = 1000.0; // Maximum distance a ray can travel before being considered invalid
float EXPOSURE = 1.0 / 8.0; // A smaller fraction makes it darker, higher makes it brighter
float EPS_ANGLE = 0.001; // Angle between each ray that gets casted from each light
float RAYS_PER_FRAME = 20; // How many rays to cast each frame
int MAX_SPINS = 3; // How many 'spins' around the origin will be executed before stopping
int DRAW_FIRST_RAY_DEPTH = 0; // Change from 0 to MAX_RAY_DEPTH, changes how many bounces a ray must have taken before it will be drawn
int RANDOM_SEED = 100; // Change to any number to make a new random scene
boolean DEBUG_DRAW = false; // Make true to see a debug view of all the rays
int CANVAS_WIDTH = 720;
int CANVAS_HEIGHT = 540;

void setup() {
  size(720, 540);
  internal_setup();

  for (int i = 0; i < 1; i++) {
    PVector pos = new PVector(200, CANVAS_HEIGHT/2 + 100);
    float radius = random(10, 100);
    Spectrum col = new Spectrum(0.5, 0.4, 0.8);
    float ior = 0.8;
    float reflectivity = 0.3;
    Circle circle = new Circle(pos, radius, col, ior, reflectivity);
    world_objects.add(circle);
  }
  for (int i = 0; i < 2; i++) {
    PVector pos = new PVector(520, CANVAS_HEIGHT/2 - 100);
    float radius = random(30, 100);
    Spectrum col = new Spectrum(0.8, 0.2, 0.3);
    float ior = 1.2;
    float reflectivity = 0.7;
    Star star = new Star(pos, radius, col, ior, reflectivity);
    world_objects.add(star);
  }

  
  lights.add(new Light(
    new PVector(CANVAS_WIDTH/2, CANVAS_HEIGHT/2),
    new Spectrum(1.0, 0.6, 0.3)
  ));

}

void draw() {
  if (DEBUG_DRAW)
    background(0);
  
  copy_image();
  
  blendMode(ADD);
  fill(25, 15, 30);
  rect(0, 0, CANVAS_WIDTH, CANVAS_WIDTH);
  blendMode(BLEND);
    
  for (int i = 0; i < RAYS_PER_FRAME; i ++) {
    if (angle > last_spin_angle + TWO_PI)
      continue;
      
    for (int n = 0; n < lights.size(); n++) {
      Light light = lights.get(n);

      ArrayList<Intersection> intersections = new ArrayList<Intersection>();
      
      Ray ray = new Ray();
      ray.origin = light.pos.copy();
      ray.dir = PVector.fromAngle(angle);
      ray.dist = 0.0;
      ray.kind = RayKind.Reflection;
      ray.depth = 0;
      ray.valid = true;
      
      ArrayList<Ray> rays = new ArrayList<Ray>();
      rays.add(ray);
      
      march(rays, intersections);
      
      Spectrum light_power = light.power.attenuate(RAY_POWER_RATIO);
      trace(intersections, 0, light_power);
    }
    
    angle += EPS_ANGLE;
  }

  if (angle >= last_spin_angle + TWO_PI) {
    angle += random(0.2);
    last_spin_angle = angle;
    spins++;
  }
  
  if (spins >= MAX_SPINS) {
    noLoop();
    save("rmpt.png");
  }
}

void keyPressed() {
  saveFrame("rmpt######.png");
}

///////////////////////////////////////////
///////////////////////////////////////////
//////////////// SHAPES ///////////////////
///////////////////////////////////////////
///////////////////////////////////////////

public class Star implements WorldObject {
  public PVector pos;
  public float radius;
  
  private Spectrum col;
  private float ior;
  private float reflectivity;
  
  public Star(PVector pos, float radius, Spectrum col, float ior, float reflectivity) {
    this.pos = pos;
    this.radius = radius;
    this.col = col;
    this.ior = ior;
    this.reflectivity = reflectivity;
  }
  
  public float get_distance(PVector pos) {
    PVector norm_pos = PVector.sub(pos, this.pos);
    float rect_d = rectangle(norm_pos, this.radius, this.radius);
    float c1_d = circle(PVector.sub(norm_pos, new PVector(-this.radius, -this.radius)), this.radius);
    float c2_d = circle(PVector.sub(norm_pos, new PVector(this.radius, -this.radius)), this.radius);
    float c3_d = circle(PVector.sub(norm_pos, new PVector(this.radius, this.radius)), this.radius);
    float c4_d = circle(PVector.sub(norm_pos, new PVector(-this.radius, this.radius)), this.radius);
    return difference(difference(difference(difference(rect_d, c1_d), c2_d), c3_d), c4_d);
  }
  
  public Spectrum evaluate_brdf(PVector light_vec, PVector normal) {
    return this.col;
  }
  
  public float get_ior() {
    return this.ior;
  }
  
  public float get_reflectivity() {
    return this.reflectivity;
  }
}

public class Circle implements WorldObject {
  public PVector pos;
  public float radius;
  
  private Spectrum col;
  private float ior;
  private float reflectivity;
  
  public Circle(PVector pos, float radius, Spectrum col, float ior, float reflectivity) {
    this.pos = pos;
    this.radius = radius;
    this.col = col;
    this.ior = ior;
    this.reflectivity = reflectivity;
  }
  
  public float get_distance(PVector pos) {
    return circle(PVector.sub(pos, this.pos), this.radius);
  }
  
  public Spectrum evaluate_brdf(PVector light_vec, PVector normal) {
    return this.col;
  }
  
  public float get_ior() {
    return this.ior;
  }
  
  public float get_reflectivity() {
    return this.reflectivity;
  }
}

//////////////////////////////////////////////////
//////////////////////////////////////////////////
//////////////// DO NOT CHANGE ///////////////////
//////////////////////////////////////////////////
//////////////////////////////////////////////////

float EPS = 0.001;
float HIT_THRESHOLD = 0.0001;
int MAX_MARCHES = 40;
float RAY_POWER_RATIO = 0.1;

float angle = 0.0;
float last_spin_angle = 0.0;
int spins = 0;
ArrayList<Light> lights;
ArrayList<WorldObject> world_objects;
float[] imgbufR;
float[] imgbufG;
float[] imgbufB;

public void internal_setup() {
  randomSeed(RANDOM_SEED);
  
  if (DEBUG_DRAW)
    RAYS_PER_FRAME = 1;

  world_objects = new ArrayList<WorldObject>();
  lights = new ArrayList<Light>();
  
  imgbufR = new float[CANVAS_WIDTH * CANVAS_HEIGHT];
  imgbufG = new float[CANVAS_WIDTH * CANVAS_HEIGHT];
  imgbufB = new float[CANVAS_WIDTH * CANVAS_HEIGHT];
}

public void copy_image() {
  for (int x = 0; x < CANVAS_WIDTH; x++) {
    for (int y = 0; y < CANVAS_HEIGHT; y++) {
      int idx = x + (y * CANVAS_WIDTH);
      float fr = imgbufR[idx] * EXPOSURE;
      float fg = imgbufG[idx] * EXPOSURE;
      float fb = imgbufB[idx] * EXPOSURE;
      fr /= 1.0 + fr;
      fg /= 1.0 + fg;
      fb /= 1.0 + fb;
      int r = min(255, (int)(fr * 255.0));
      int g = min(255, (int)(fg * 255.0));
      int b = min(255, (int)(fb * 255.0));
      set(x, y, color(r, g, b));
    }
  }
}

float intersection(float a, float b) {
    return max(b, b);
}

float union(float a, float b) {
  return min(a, b);
}

float difference(float a, float b) {
  return max(a, -b);
}

public float rectangle(PVector pos, float w, float h) {
  PVector component_dist = new PVector(abs(pos.x) - w, abs(pos.y) - h);
  float outside_dist = new PVector(max(component_dist.x, 0), max(component_dist.y, 0)).mag();
  float inside_dist = min(max(component_dist.x, component_dist.y), 0);
  return outside_dist + inside_dist;
}

public float circle(PVector pos, float r) {
  return pos.mag() - r;
}

public float fresnel(float reflectivity, PVector l, PVector n) {
  float ndotl = max(0, PVector.dot(l, n));
  return reflectivity + (1.0 - reflectivity) * pow(1 - ndotl, 5);
}

public void trace(ArrayList<Intersection> intersections, int curr_intersection_id, Spectrum light_power) {
  Intersection intersection = intersections.get(curr_intersection_id);
  if (intersection.ray.valid) {
    PVector origin = intersection.ray.origin;
    PVector end = PVector.add(intersection.ray.origin, PVector.mult(intersection.ray.dir, intersection.ray.dist));
    
    if (intersection.ray.depth >= DRAW_FIRST_RAY_DEPTH)
      drawLine(origin.x, origin.y, end.x, end.y, light_power);
    
    if (DEBUG_DRAW) {
      if (intersection.ray.kind == RayKind.Reflection) {
        stroke(0, 255, 0);
      } else {
        stroke(0, 0, 255);
      }
      strokeWeight(1);
      line(origin.x, origin.y, end.x, end.y);
    }
    
    PVector light_vec = PVector.mult(intersection.ray.dir, -1.0);
    WorldObject obj = world_objects.get(intersection.obj_id.id);
    Spectrum brdf = obj.evaluate_brdf(light_vec, intersection.normal);
    
    float fres = fresnel(obj.get_reflectivity(), light_vec, intersection.normal);
    
    if (intersection.ray.kind == RayKind.Refraction) {
      fres = 1.0 - fres;
    }
    
    if (intersection.ray.depth < MAX_RAY_DEPTH) {
      trace(intersections, curr_intersection_id * 2 + 1, light_power.attenuate(fres).attenuate(brdf));
      trace(intersections, curr_intersection_id * 2 + 2, light_power.attenuate(fres).attenuate(brdf));
    }
  }
}

public class WorldObjectId {
  public int id;
  
  public WorldObjectId(int id) {
    this.id = id;
  }
}

public class Intersection {
  public Ray ray;
  public PVector normal;
  public WorldObjectId obj_id;
}

public enum RayKind {
  Reflection,
  Refraction
}

public class Ray {
  public PVector origin;
  public PVector dir;
  public float dist;
  public int depth;
  public RayKind kind;
  public boolean valid;
}

public class Light {
  public PVector pos;
  public Spectrum power;
  
  public Light(PVector pos, Spectrum power) {
    this.pos = pos;
    this.power = power;
  }
}

public class Spectrum {
  public float r;
  public float g;
  public float b;
  
  public Spectrum(float r, float g, float b) {
    this.r = r;
    this.g = g;
    this.b = b;
  }
  
  public int r() {
    return (int)(this.r * 255.0);
  }
  public int g() {
    return (int)(this.g * 255.0);
  }
  public int b() {
    return (int)(this.b * 255.0);
  }
  
  public Spectrum reinhard() {
    Spectrum new_spect = new Spectrum(this.r, this.g, this.b);
    new_spect.r /= (1.0 + this.r);
    new_spect.g /= (1.0 + this.g);
    new_spect.b /= (1.0 + this.b);
    return new_spect;
  }
  
  public Spectrum attenuate(float scalar) {
    Spectrum new_spect = new Spectrum(this.r, this.g, this.b);
    new_spect.r *= scalar;
    new_spect.g *= scalar;
    new_spect.b *= scalar;
    return new_spect;
  }
  
  public Spectrum attenuate(Spectrum spect) {
    Spectrum new_spect = new Spectrum(this.r, this.g, this.b);
    new_spect.r *= spect.r;
    new_spect.g *= spect.g;
    new_spect.b *= spect.b;
    return new_spect;
  }
}

public interface SDF {
  public float get_distance(PVector pos);
}

public interface BRDF {
  public Spectrum evaluate_brdf(PVector incident, PVector normal);
}

public interface Fresnel {
  public float get_reflectivity();
}

public interface IOR {
  public float get_ior();
}

public interface WorldObject extends SDF, BRDF, Fresnel, IOR {
}

PVector reflect(PVector incident, PVector normal) {
  return PVector.sub(incident, PVector.mult(normal, 2 * PVector.dot(incident, normal)));
}

PVector refract(PVector incident, PVector normal, float ior) {
  float ndi = PVector.dot(normal, incident);
  float k = 1.0 - ior * ior * (1.0 - ndi * ndi);
  if (k < 0.0)
    return new PVector(0.0, 0.0);
  else
    return PVector.sub(PVector.mult(incident, ior), PVector.mult(normal, ior * ndi + sqrt(k)));
}

float dist_to_scene(PVector p, WorldObjectId closest) {
  float min_distance = MAX_FLOAT;
  for (int i = 0; i < world_objects.size(); i++) {
    WorldObject obj = world_objects.get(i);
    float distance = obj.get_distance(p);
    if (distance < min_distance) {
      min_distance = distance;
      closest.id = i;
    }
  }
  return min_distance;
}

PVector estimate_normal(PVector p){
  WorldObjectId id = new WorldObjectId(0);
  float x_p = dist_to_scene(new PVector(p.x+EPS,p.y), id);
  float x_m = dist_to_scene(new PVector(p.x-EPS,p.y), id);
  float y_p = dist_to_scene(new PVector(p.x,p.y+EPS), id);
  float y_m = dist_to_scene(new PVector(p.x,p.y-EPS), id);
  float x_diff = x_p - x_m;
  float y_diff = y_p - y_m;
  PVector vec = new PVector(x_diff, y_diff);
  vec.normalize();
  return vec;
}

void march(
  ArrayList<Ray> rays,
  ArrayList<Intersection> intersections
) {
  while (rays.size() > 0) {
    Ray ray = rays.remove(0);
    
    WorldObjectId closest = new WorldObjectId(0);
    float min_distance = MAX_FLOAT;
    float prev_min_distance = MAX_FLOAT;
    int marches = 0;
    
    PVector pos = new PVector(0, 0);
    
    if (ray.valid) {
      while (abs(min_distance) > HIT_THRESHOLD && marches < MAX_MARCHES && ray.dist < MAX_RAY_DIST) {
        pos = PVector.add(ray.origin, PVector.mult(ray.dir, ray.dist));
        prev_min_distance = min_distance;
        min_distance = dist_to_scene(pos.copy(), closest);
        ray.dist += abs(min_distance);
        marches++;
      }
    }
    
    Intersection intersection = new Intersection();
    intersection.ray = ray;
    if (ray.valid)
      intersection.normal = estimate_normal(pos);
    if (prev_min_distance < 0)
      intersection.normal.mult(-1.0);
    if (DEBUG_DRAW && ray.valid) {
      stroke(255, 0, 0);
      strokeWeight(1);
      line(pos.x, pos.y, pos.x + intersection.normal.x * 10, pos.y + intersection.normal.y * 10);
    }
    intersection.obj_id = closest;
  
    if (ray.depth < MAX_RAY_DEPTH) {
      if (ray.valid && ray.dist < MAX_RAY_DIST) {
        Ray reflection_ray = new Ray();
        reflection_ray.dir = reflect(ray.dir, intersection.normal);
        reflection_ray.origin = PVector.add(pos, PVector.mult(reflection_ray.dir, 0.01));
        reflection_ray.dist = 0.0;
        reflection_ray.kind = RayKind.Reflection;
        reflection_ray.depth = ray.depth + 1;
        reflection_ray.valid = true;
        rays.add(reflection_ray);
      
        WorldObject obj = world_objects.get(closest.id);
        Ray refraction_ray = new Ray();
        refraction_ray.dir = refract(ray.dir, intersection.normal, obj.get_ior());
        refraction_ray.origin = PVector.add(pos, PVector.mult(refraction_ray.dir, 0.01));
        refraction_ray.dist = 0.0;
        refraction_ray.kind = RayKind.Refraction;
        refraction_ray.depth = ray.depth + 1;
        refraction_ray.valid = true;
        rays.add(refraction_ray);
      } else {
        Ray reflection_ray = new Ray();
        reflection_ray.valid = false;
        reflection_ray.depth = ray.depth + 1;
        rays.add(reflection_ray);
        Ray refraction_ray = new Ray();
        refraction_ray.valid = false;
        refraction_ray.depth = ray.depth + 1;
        rays.add(refraction_ray);
      }
    }
    
    intersections.add(intersection);
  }
}

void plot(double xd, double yd, double a, Spectrum c) {
  int x = round((float)xd);
  int y = round((float)yd);
  if (x >= CANVAS_WIDTH || x < 0 || y >= CANVAS_HEIGHT || y < 0) {
    return;
  }
  int idx = x + (CANVAS_WIDTH * y);
  Spectrum mapped = c.attenuate((float)a);
  
  imgbufR[idx] += mapped.r;
  imgbufG[idx] += mapped.g;
  imgbufB[idx] += mapped.b;
}

int ipart(double x) {
    return (int)x;
}
 
double fpart(double x) {
    return x - Math.floor(x);
}
 
double rfpart(double x) {
    return 1.0 - fpart(x);
}
 
void drawLine(float x0, float y0, float x1, float y1, Spectrum s) {
    boolean steep = Math.abs(y1 - y0) > Math.abs(x1 - x0);
    if (steep) {
      float t = x0;
      x0 = y0;
      y0 = t;
      t = x1;
      x1 = y1;
      y1 = t;
    }
 
    if (x0 > x1) {
      float t = x0;
      x0 = x1;
      x1 = t;
      t = y0;
      y0 = y1;
      y1 = t;
    }
 
    double dx = x1 - x0;
    double dy = y1 - y0;
    double gradient = dy / dx;
    
    double angle = Math.atan2(dy, dx);
    double m = 0.5 * Math.abs(Math.sin(angle * 2)) + 0.5;

    // handle first endpoint
    double xend = Math.round(x0);
    double yend = y0 + gradient * (xend - x0);
    double xgap = rfpart(x0 + 0.5);
    double xpxl1 = xend; // this will be used in the main loop
    double ypxl1 = ipart(yend);
 
    if (steep) {
        plot(ypxl1, xpxl1, rfpart(yend) * xgap * m, s);
        plot(ypxl1 + 1, xpxl1, fpart(yend) * xgap * m, s);
    } else {
        plot(xpxl1, ypxl1, rfpart(yend) * xgap * m, s);
        plot(xpxl1, ypxl1 + 1, fpart(yend) * xgap * m, s);
    }
 
    // first y-intersection for the main loop
    double intery = yend + gradient;
 
    // handle second endpoint
    xend = Math.round(x1);
    yend = y1 + gradient * (xend - x1);
    xgap = fpart(x1 + 0.5);
    double xpxl2 = xend; // this will be used in the main loop
    double ypxl2 = ipart(yend);
 
    if (steep) {
        plot(ypxl2, xpxl2, rfpart(yend) * xgap * m, s);
        plot(ypxl2 + 1, xpxl2, fpart(yend) * xgap * m, s);
    } else {
        plot(xpxl2, ypxl2, rfpart(yend) * xgap * m, s);
        plot(xpxl2, ypxl2 + 1, fpart(yend) * xgap * m, s);
    }
 
    // main loop
    for (double x = xpxl1 + 1; x <= xpxl2 - 1; x++) {
        if (steep) {
            plot(ipart(intery), x, rfpart(intery) * m, s);
            plot(ipart(intery) + 1, x, fpart(intery) * m, s);
        } else {
            plot(x, ipart(intery), rfpart(intery) * m, s);
            plot(x, ipart(intery) + 1, fpart(intery) * m, s);
        }
        intery = intery + gradient;
    }
}
