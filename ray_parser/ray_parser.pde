///////////////////////////////////////////////////////////////////////
//
// Command Line Interface (CLI) Parser  
//
///////////////////////////////////////////////////////////////////////
String gCurrentFile = new String("rect_test.cli"); // A global variable for holding current active file name.

///////////////////////////////////////////////////////////////////////
//
// Press key 1 to 9 and 0 to run different test cases.
//
///////////////////////////////////////////////////////////////////////
void keyPressed() {
  clearAll();
  switch(key) {
    case '1':  gCurrentFile = new String("t0.cli"); interpreter(); break;
    case '2':  gCurrentFile = new String("t1.cli"); interpreter(); break;
    case '3':  gCurrentFile = new String("t2.cli"); interpreter(); break;
    case '4':  gCurrentFile = new String("t3.cli"); interpreter(); break;
    case '5':  gCurrentFile = new String("c0.cli"); interpreter(); break;
    case '6':  gCurrentFile = new String("c1.cli"); interpreter(); break;
    case '7':  gCurrentFile = new String("c2.cli"); interpreter(); break;
    case '8':  gCurrentFile = new String("c3.cli"); interpreter(); break;
    case '9':  gCurrentFile = new String("c4.cli"); interpreter(); break;
    case '0':  gCurrentFile = new String("c5.cli"); interpreter(); break;
  }
}

///////////////////////////////////////////////////////////////////////
//
//  Parser core. It parses the CLI file and processes it based on each 
//  token. Only "color", "rect", and "write" tokens are implemented. 
//  You should start from here and add more functionalities for your
//  ray tracer.
//
//  Note: Function "splitToken()" is only available in processing 1.25 
//       or higher.
//
///////////////////////////////////////////////////////////////////////

final float SCRDIM = 300;
final int MAXDEPTH = 20;
float fov = 45.;
color bg;
Surface curSurface;
Triangle curTriangle;
ArrayList<Light> lights;
ArrayList<Sphere> spheres;
ArrayList<Triangle> triangles;

void clearAll()
{
  fill(0, 0, 0);
  rect(0, 0, SCRDIM, SCRDIM);
  
  setAll();
}

void setAll()
{
    lights = new ArrayList<Light>();
  spheres = new ArrayList<Sphere>();
  triangles = new ArrayList<Triangle>();
}

void interpreter() {
  
  String str[] = loadStrings(gCurrentFile);
  if (str == null) println("Error! Failed to read the file.");
  for (int i=0; i<str.length; i++) {
    
    String[] token = splitTokens(str[i], " "); // Get a line and parse tokens.
    if (token.length == 0) continue; // Skip blank line.
    
    if (token[0].equals("fov")) {
      fov =float(token[1]);
    }
    else if (token[0].equals("background")) {
      float bgr =float(token[1]);
      float bgg =float(token[2]);
      float bgb =float(token[3]);
      bg = color(bgr, bgg, bgb);
    }
    else if (token[0].equals("light")) {
      float x =float(token[1]);
      float y =float(token[2]);
      float z =float(token[3]);
      float r =float(token[4]);
      float g =float(token[5]);
      float b =float(token[6]);
      lights.add(new Light(x, y, z, r, g, b));
    }
    else if (token[0].equals("surface")) {
      float r1 =float(token[1]);
      float g1 =float(token[2]);
      float b1 =float(token[3]);
      float r2 =float(token[4]);
      float g2 =float(token[5]);
      float b2 =float(token[6]);
      float r3 =float(token[7]);
      float g3 =float(token[8]);
      float b3 =float(token[9]);
      float p =float(token[10]);
      float k =float(token[11]);
      curSurface = new Surface(r1, g1, b1, r2, g2, b2, r3, g3, b3, p, k);
    }    
    else if (token[0].equals("sphere")) {
      float r =float(token[1]);
      float x =float(token[2]);
      float y =float(token[3]);
      float z =float(token[4]);
      Sphere s = new Sphere(x, y, z, r);
      s.surface = curSurface;
      spheres.add(s);
    }
    else if (token[0].equals("begin")) {
      curTriangle = new Triangle();
      curTriangle.surface = curSurface;
    }
    else if (token[0].equals("vertex")) {
      float x =float(token[1]);
      float y =float(token[2]);
      float z =float(token[3]);
      curTriangle.addPoint(x, y, z);
    }
    else if (token[0].equals("end")) {
      triangles.add(curTriangle);
    }
    else if (token[0].equals("color")) {
      float r =float(token[1]);
      float g =float(token[2]);
      float b =float(token[3]);
      fill(r, g, b);
    }
    else if (token[0].equals("rect")) {
      float x0 = float(token[1]);
      float y0 = float(token[2]);
      float x1 = float(token[3]);
      float y1 = float(token[4]);
      rect(x0, height-y1, x1-x0, y1-y0);
    }
    else if (token[0].equals("write")) {
      
      float k = tan(radians(fov / 2));
      
      fill(100, 0, 0);
      if (triangles.size() + spheres.size() > 0)
      {
        for (float a = 0; a < SCRDIM; a++)
        {
          for (float b = 0; b < SCRDIM; b++)
          {
            
            
            float ta = (a - SCRDIM / 2) * (2 * k / SCRDIM);
            float tb = (SCRDIM / 2 - b) * (2 * k / SCRDIM);
            float[][] morig = {{0,0,0}};
            float[][] mray = {{ta, tb, -1}};
            gtMatrix orig = new gtMatrix(morig);
            gtMatrix ray = new gtMatrix(mray);
            
            fill(rayCast(orig, ray, 0, null));
            rect(a, b, 1, 1);
          }
        }
      }  
      
      save(token[1]);  
    }
  }
}

// my lib

// helper fx

color rayCast(gtMatrix orig, gtMatrix ray, int depth, GeoPrimitive avoid)
{
  Collision coll = objHit(orig, ray, avoid);
  
  
  if (coll  == null)
  {
    return bg;
  }
  else
  {
    Surface s = coll.gp.surface;
    color diffuse = s.diffuse;
    color ambient = s.ambient;
    color specular = s.specular;
    
    float rsum = red(ambient);
    float gsum = green(ambient);
    float bsum = blue(ambient);
    gtMatrix eyeray = orig.subtract(coll.loc).normalize();
    for(Light l : lights)
    { 
      
      
      gtMatrix lray = l.loc.subtract(coll.loc).normalize();
      gtMatrix reflray = coll.norm.multiply(2).multiply(coll.norm.dotProduct(lray)).subtract(lray);      
      Collision coll2 = objHit(l.loc, lray.multiply(-1.0), null);

      if ((coll2 != null) && (coll2.gp == coll.gp))
      {
        float max1 = max(0, coll.norm.dotProduct(lray));
        float max2 = max(0, eyeray.dotProduct(reflray));
  
        rsum += red(l.col) * (red(diffuse) * max1 + red(specular) * pow(max2, s.phong));
        gsum += green(l.col) * (green(diffuse) * max1 + green(specular) * pow(max2, s.phong));
        bsum += blue(l.col) * (blue(diffuse) * max1 + blue(specular) * pow(max2, s.phong));
        
      }
    }
    float refl = s.refl;
    
    gtMatrix reflray = coll.norm.multiply(2).multiply(coll.norm.dotProduct(eyeray)).subtract(eyeray);      
    
    if ((refl > 0) && (depth < MAXDEPTH))
    {
      color color2 = rayCast(coll.loc, reflray, depth + 1, coll.gp);
      float rsum2 = red(color2);
      float gsum2 = green(color2);
      float bsum2 = blue(color2);
      return color(rsum + rsum2 * refl, gsum + gsum2 * refl, bsum + bsum2 * refl); 
    }
    return color(rsum, gsum, bsum);
  }
   
}

Collision objHit(gtMatrix orig, gtMatrix ray, GeoPrimitive avoid)
{
  Collision c1 = sphereHit(orig, ray, avoid);
  Collision c2 = triangleHit(orig, ray, avoid);
  if (c1 == null)
    return c2;
  if (c2 == null)
    return c1;
  if (c1.dist < c2.dist)
    return c1;
  return c2;
}

Collision triangleHit(gtMatrix orig, gtMatrix ray, GeoPrimitive avoid)
{
  float tdist = -1;
  Triangle closest = null;
  for (Triangle tr : triangles)
  {
    if ( (avoid != null) && (tr == avoid) )
      continue;
    gtMatrix qo = orig.subtract(tr.center);
    gtMatrix norm = cross3(tr.v1.subtract(tr.v2), tr.v1.subtract(tr.v3));
    float t = - qo.dotProduct(norm) / ray.dotProduct(norm);
    gtMatrix loc = ray.multiply(t).add(orig);
    if (sameSide(tr.v1, tr.v2, tr.v3, loc) && sameSide(tr.v2, tr.v3, tr.v1, loc) && sameSide(tr.v3, tr.v1, tr.v2, loc))
    {
      if (t > 0 && (tdist > t || (tdist == -1)))
      {
        closest = tr;
        tdist = t;
      }
    }
  }
  if (tdist == -1)
  {
    return null;
  }
  gtMatrix loc = ray.multiply(tdist).add(orig);
  gtMatrix norm = cross3(closest.v2.subtract(closest.v1), closest.v1.subtract(closest.v3));
  gtMatrix disp = loc.subtract(orig);
  return new Collision(closest, loc, ray, norm, disp.magnitude());
}
  
    
boolean sameSide(gtMatrix a, gtMatrix b, gtMatrix c, gtMatrix p)
{
  return (cross3(b.subtract(a), c.subtract(a)).dotProduct(cross3(b.subtract(a), p.subtract(a))) > 0);
}

Collision sphereHit(gtMatrix orig, gtMatrix ray, GeoPrimitive avoid)
{
  Sphere closest = null;
  float tdist = -1;
  
  for (int i = 0; i < spheres.size(); i++)
  {    
    Sphere s = spheres.get(i);
    if ( (avoid != null) && (s == avoid) )
      continue;
    gtMatrix center = s.center;
    gtMatrix co = orig.subtract(center);
        
    float a = ray.dotProduct(ray);
    float b = ray.multiply(2).dotProduct(co);
    float c = co.dotProduct(co) - s.radius * s.radius;
    
    float disc = b * b - 4 * a * c;
    if (disc < 0)
    {
      continue;
    }
    else
    {
      float t = min((-b - sqrt(disc)) / (2 * a), (-b + sqrt(disc)) / (2 * a)) ;
      if (t > 0 && (tdist > t || (tdist == -1)))
      {
        closest = s;
        tdist = t;
      }
    }
  }
  if (tdist == -1)
  {
    return null;
  }
  gtMatrix loc = ray.multiply(tdist).add(orig);
  gtMatrix norm = loc.subtract(closest.center);
  gtMatrix disp = loc.subtract(orig);
  return new Collision(closest, loc, ray, norm, disp.magnitude());
}


// data classes

class Collision
{
  public GeoPrimitive gp = null;
  public gtMatrix loc, inc, norm;
  public float dist;
  public Collision(GeoPrimitive cgp, gtMatrix cloc, gtMatrix cinc, gtMatrix cnorm, float cdist)
  {
    gp = cgp;
    loc = cloc;
    inc = cinc.normalize();
    norm = cnorm.normalize();
    dist = cdist;
  }
}

abstract class GeoPrimitive
{
  public Surface surface;
}

class Sphere extends GeoPrimitive
{
  public float radius;
  public gtMatrix center;
  public Sphere(float sx, float sy, float sz, float sr)
  {
    float[][] c = {{sx, sy, sz}};
    center = new gtMatrix(c);
    radius = sr;
  }
}

class Triangle extends GeoPrimitive
{
  gtMatrix v1, v2, v3, center;
  int curPoint = 0;
  public void addPoint(float tx, float ty, float tz)
  {
    float[][] v = {{tx, ty, tz}};
    curPoint++;
    switch (curPoint)
    {
      case 1:
        v1 = new gtMatrix(v);
        break;
      case 2:
        v2 = new gtMatrix(v);
        break;
      case 3:
        v3 = new gtMatrix(v);
        center = v1.add(v2).add(v3).multiply(1.0 / 3.0);
        break;
    }
  }
}

class Light
{
  public gtMatrix loc;
  public color col;
  public Light(float lx, float ly, float lz, float cr, float cg, float cb)
  {
    float[][] m2 = {{lx, ly, lz}};
    loc = new gtMatrix(m2);
    col = color(cr, cg, cb);
  }
}

class Surface
{
  public color diffuse, ambient, specular;
  public float phong, refl;
  
  public Surface(float cdr, float cdg, float cdb, float car, float cag, float cab, float csr, float csg, float csb, float cp, float ck)
  {
    diffuse = color(cdr, cdg, cdb);
    ambient = color(car, cag, cab);
    specular = color(csr, csg, csb);
    phong = cp;
    refl = ck;
  }
}

class gtMatrix {
  float[][] m;
  public int i, j;
  gtMatrix(int mi, int mj) 
  {
    i = mi; j = mj;
    m = new float[i][j];
  }
  gtMatrix(float[][] m2) 
  {
    setMatrix(m2);
  }
  void identity()
  {
    for(int mi = 0; mi < i; mi++)
    {
      for (int mj = 0; mj < j; mj++)
      {
        if (mi == mj)
          setVal(mi, mj, 1);
        else
          setVal(mi, mj, 0);
      }
    }
  }
  void setVal(int mi, int mj, float val)
  {
    m[mi][mj] = val;
  }
  float getVal(int mi, int mj)
  {
    return m[mi][mj];
  }
  void setMatrix(float[][] m2)
  {
    m = m2;
    i = m2.length;
    j = m2[0].length;
  }
  float[][] getMatrix()
  {
    return m;
  }
  gtMatrix add(gtMatrix m2)
  {
    return gtMatAdd(this, m2);
  }
  gtMatrix subtract(gtMatrix m2)
  {
    return gtMatAdd(this, m2.multiply(-1.));
  }
  gtMatrix multiply(gtMatrix m2)
  {
    return gtMatMult(this, m2);
  }
  gtMatrix multiply(float f)
  {
    return gtMatConstMult(this, f);
  }
  gtMatrix transpose()
  {
    return gtMatTranspose(this);
  }
  gtMatrix normalize()
  {
    return gtMatNorm(this);
  }
  float magnitude()
  {
    return gtMatMagn(this);
  }
  float dotProduct(gtMatrix m2)
  {
    return this.multiply(gtMatTranspose(m2)).getVal(0,0);
  }
  
  String toString()
  {
    String s = "[";
    for (int si = 0; si < i; si++)
    {
      s += "[";
      for (int sj = 0; sj < j; sj++)
      {
        s += getVal(si, sj);
        if (sj + 1 < j)
          s += " ";
      }
      if (si + 1 < i)
        s += "]\n ";
      else
        s += "]]";
    }
    return s;
  } 
}

gtMatrix gtMatAdd(gtMatrix a, gtMatrix b)
{
  gtMatrix res = new gtMatrix(a.i, a.j);
  for (int i = 0; i < a.i; i++)
  {
    for (int j = 0; j < a.j; j++)
    {
      res.setVal(i, j, a.getVal(i, j) + b.getVal(i,j));
    }
  }
  return res;
}

gtMatrix gtMatConstMult(gtMatrix a, float f)
{
  gtMatrix res = new gtMatrix(a.i, a.j);
  for (int i = 0; i < a.i; i++)
  {
    for (int j = 0; j < a.j; j++)
    {
      res.setVal(i, j, a.getVal(i, j) * f);
    }
  }
  return res;
}

gtMatrix gtMatMult(gtMatrix a, gtMatrix b)
{
  gtMatrix res = new gtMatrix(a.i, b.j);
  for (int i = 0; i < a.i; i++)
  {
    for (int j = 0; j < b.j; j++)
    {
      float sum = 0;
      for (int k = 0; k < a.j; k++)
      {
        sum += a.getVal(i, k) * b.getVal(k, j);
      }
      res.setVal(i, j, sum);
    }
  }
  return res;
}


gtMatrix gtMatTranspose(gtMatrix a)
{
  gtMatrix res = new gtMatrix(a.j, a.i);
  for (int i = 0; i < a.i; i++)
  {
    for (int j = 0; j < a.j; j++)
    {
      res.setVal(j, i, a.getVal(i, j));
    }
  }
  return res;
}

gtMatrix gtMatNorm(gtMatrix a)
{
  gtMatrix res = new gtMatrix(a.i, a.j);
  res.setMatrix(a.m);
  res = res.multiply(1.0 / gtMatMagn(a));
  return res;
}

float gtMatMagn(gtMatrix a)
{
  float sum = 0;
  gtMatrix res = new gtMatrix(a.i, a.j);
  for (int i = 0; i < a.i; i++)
  {
    for (int j = 0; j < a.j; j++)
    {
      float val = a.getVal(i, j);
      sum += val * val;
      res.setVal(i, j, val);
    }
  }
  sum = sqrt(sum);
  return sum;
}

gtMatrix cross3(gtMatrix m1, gtMatrix m2)
{
  float u1 = m1.getVal(0,0);
  float u2 = m1.getVal(0,1);
  float u3 = m1.getVal(0,2);
  float v1 = m2.getVal(0,0);
  float v2 = m2.getVal(0,1);
  float v3 = m2.getVal(0,2);
  float[][] m3 = {{u2 * v3 - u3 * v2, u3 * v1 - u1 * v3, u1 * v2 - u2 * v1}};
  return new gtMatrix(m3);
}

///////////////////////////////////////////////////////////////////////
//
// Some initializations for the scene.
//
///////////////////////////////////////////////////////////////////////
void setup() {
  size(300, 300);  
  noStroke();
  colorMode(RGB, 1.0);
  background(0, 0, 0);
  setAll();
  interpreter();
}

///////////////////////////////////////////////////////////////////////
//
// Draw frames.  Should leave this empty.
//
///////////////////////////////////////////////////////////////////////
void draw() {
}

