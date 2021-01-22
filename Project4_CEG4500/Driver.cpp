// Driver.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include <iostream>
#include <GL/glut.h>
#include <fstream>
#include <cstdio>

#include "ply.h"
#include <vector>
using namespace std;

class Vec {
public:
    float x, y, z;

    Vec() {
        x = y = z = 0;
    }

    Vec(const Vec& v2) {
        this->x = v2.x;
        this->y = v2.y;
        this->z = v2.z;
    }

    Vec(float x, float y, float z) {
        this->x = x;
        this->y = y;
        this->z = z;
    }

    Vec( Vertex v) {
        this->x = v.x;
        this->y = v.y;
        this->z = v.z;
    }

    //add
    Vec operator + (Vec const& obj) {
        Vec res(*this);
        res.x += obj.x;
        res.y += obj.y;
        res.z += obj.z;
        return res;
    }

    //subtract
    Vec operator - (Vec const& obj) {
        Vec res(*this);
        res.x -= obj.x;
        res.y -= obj.y;
        res.z -= obj.z;
        return res;
    }

    //dot product
    float operator * (Vec const& obj) {
        return this->x * obj.x + this->y * obj.y + this->z * obj.z;
    }

    //cross product
    Vec operator / (Vec const& obj) {
        return Vec(
            this->y * obj.z - this->z * obj.y,
            -(this->x * obj.z - this->z * obj.x),
            this->x * obj.y - this->y * obj.x
        );
    }

    //scalar multiply
    Vec operator * (float const& obj) {
        Vec res(*this);
        res.x *= obj;
        res.y *= obj;
        res.z *= obj;
        return res;
    }
    //scalar divide
    Vec operator / (float const& obj) {
        Vec res(*this);
        res.x /= obj;
        res.y /= obj;
        res.z /= obj;
        return res;
    }


    float mag() {
        return sqrtf(powf(x, 2) + powf(y, 2) + powf(z, 2));
    }

    void normalize() {
        float magnitude = mag();

        this->x /= magnitude;
        this->y /= magnitude;
        this->z /= magnitude;

    }
};

ostream& operator<<(ostream& os, const Vec& v)
{
    os << "(" << v.x << ", " << v.y << ", " << v.z << ")";
    return os;
}

std::vector<Vec> surface_normals;
std::vector<float> surface_area;

float fuzz = 1.0e-5;

struct camera {
    Vec eye, vd, up;
    float w, h, dist;
};

struct light {
    Vec i_i, i_a;
    float k_d, k_s, k_a;
};

camera cam;
light li;

int num_elems, num_verts, num_faces;
Vertex** vlist;
Face** flist;

// global for the longest line across object
float internal_diagonal;
// global for the corners of the bounding box
Vertex min_corner, max_corner;
// width and height of window, should automatically ajust other parts of program
const int WIDTH = 100, HEIGHT = 100;
unsigned int pixels[HEIGHT][WIDTH][3];
// ply file to read in
const char* FILE_NAME = "./icosahedron.ply";


void save_ppm(const char* filename) {
    ofstream ofs(filename, ios_base::out | ios_base::binary);
    ofs << "P6" << endl << WIDTH << ' ' << HEIGHT << endl << "255" << endl;

    for (int j = HEIGHT - 1; j >= 0; j--)
        for (int i = 0; i < WIDTH; i++)
            ofs << (char)((pixels[j][i][0] / (256*256*256)) % 256) << (char)((pixels[j][i][1] / (256 * 256 * 256)) % 256) << (char)((pixels[j][i][2] / (256 * 256 * 256)) % 256);       // red, green, blue

    ofs.close();
}


// deleted everything that printed stuff, I wanted to speed up my process.
// I also added something to figure out the bounding box around the object as I'm reading it in, rather than searching through after its all read in.
void read_test(char* filename)
{
    int i, j, k;
    PlyFile* ply;
    int nelems;
    char** elist;
    int file_type;
    float version;
    int nprops;
    PlyProperty** plist;
    char* elem_name;
    int num_comments;
    char** comments;
    int num_obj_info;
    char** obj_info;

    /* open a PLY file for reading */
    ply = ply_open_for_reading(filename, &nelems, &elist, &file_type, &version);

    /* go through each kind of element that we learned is in the file */
    /* and read them */
    for (i = 0; i < nelems; i++) {
        /* get the description of the first element */
        elem_name = elist[i];
        plist = ply_get_element_description(ply, elem_name, &num_elems, &nprops);

        /* if we're on vertex elements, read them in */
        if (equal_strings((char*)"vertex", elem_name)) {

            /* create a vertex list to hold all the vertices */
            vlist = (Vertex**)malloc(sizeof(Vertex*) * num_elems);
            num_verts = num_elems;

            /* set up for getting vertex elements */

            ply_get_property(ply, elem_name, &vert_props[0]);
            ply_get_property(ply, elem_name, &vert_props[1]);
            ply_get_property(ply, elem_name, &vert_props[2]);

            /* grab all the vertex elements */
            for (j = 0; j < num_elems; j++) {

                /* grab and element from the file */
                vlist[j] = (Vertex*)malloc(sizeof(Vertex));
                ply_get_element(ply, (void*)vlist[j]);

                /* find bounding box for ply file */
                if (j == 0) {
                    memcpy(&min_corner, vlist[j], sizeof(Vertex));
                    memcpy(&max_corner, vlist[j], sizeof(Vertex));
                }
                else {
                    if (min_corner.x > vlist[j]->x) min_corner.x = vlist[j]->x;
                    if (min_corner.y > vlist[j]->y) min_corner.y = vlist[j]->y;
                    if (min_corner.z > vlist[j]->z) min_corner.z = vlist[j]->z;
                    if (max_corner.x < vlist[j]->x) max_corner.x = vlist[j]->x;
                    if (max_corner.y < vlist[j]->y) max_corner.y = vlist[j]->y;
                    if (max_corner.z < vlist[j]->z) max_corner.z = vlist[j]->z;
                }

            }
        }

        /* if we're on face elements, read them in */
        if (equal_strings((char*)"face", elem_name)) {

            /* create a list to hold all the face elements */
            flist = (Face**)malloc(sizeof(Face*) * num_elems);
            num_faces = num_elems;

            /* set up for getting face elements */

            ply_get_property(ply, elem_name, &face_props[0]);
            ply_get_property(ply, elem_name, &face_props[1]);

            /* grab all the face elements */
            for (j = 0; j < num_elems; j++) {

                /* grab and element from the file */
                flist[j] = (Face*)malloc(sizeof(Face));
                ply_get_element(ply, (void*)flist[j]);
            }
        }
    }

    /* grab and print out the comments in the file */
    comments = ply_get_comments(ply, &num_comments);

    /* grab and print out the object information */
    obj_info = ply_get_obj_info(ply, &num_obj_info);

    /* close the PLY file */
    ply_close(ply);
}

void calc_normals() {
    // ensure normal vector is set to zero
    for (int i = 0; i < num_verts; i++) {
        vlist[i]->nx = 0;
        vlist[i]->ny = 0;
        vlist[i]->nz = 0;
    }

    for (int i = 0; i < num_faces; i++) {
        // newells method for calculating surface vectors
        float nx = 0, ny = 0, nz = 0;
        for (int j = 0; j < flist[i]->nverts; j++) {
            Vertex* current = vlist[flist[i]->verts[j]];
            Vertex* next = vlist[flist[i]->verts[(j + 1) % flist[i]->nverts]];

            nx += (current->y - next->y) * (current->z + next->z);
            ny += (current->z - next->z) * (current->x + next->x);
            nz += (current->x - next->x) * (current->y + next->y);
        }

        float length = sqrtf(powf(nx, 2) + powf(ny, 2) + powf(nz, 2));

        nx /= length;
        ny /= length;
        nz /= length;

        surface_normals.push_back(Vec(nx, ny, nz));

        // add surface normal vector to vertex normal vector
        for (int j = 0; j < flist[i]->nverts; j++) {
            vlist[flist[i]->verts[j]]->nx += nx;
            vlist[flist[i]->verts[j]]->ny += ny;
            vlist[flist[i]->verts[j]]->nz += nz;
        }
    }

    // normalize vertex normal vectors
    for (int i = 0; i < num_verts; i++) {
        float length = sqrtf(powf(vlist[i]->nx, 2) + powf(vlist[i]->ny, 2) + powf(vlist[i]->nz, 2));

        vlist[i]->nx /= length;
        vlist[i]->ny /= length;
        vlist[i]->nz /= length;
    }
}

void calc_face_area() {
    
    for (int i = 0; i < num_faces; i++) {
        Vec center;
        for (int j = 0; j < flist[i]->nverts; j++) {
            center = center + Vec(*vlist[flist[i]->verts[j]]);
        }
        center = center / flist[i]->nverts;
        float area = ((Vec(*vlist[flist[i]->verts[flist[i]->nverts - 1]]) - center) / (Vec(*vlist[flist[i]->verts[0]]) - center)).mag() / 2;
        for (int j = 0; j < flist[i]->nverts - 1; j++) {
            area += ((Vec(*vlist[flist[i]->verts[j]]) - center) / (Vec(*vlist[flist[i]->verts[j + 1]]) - center)).mag() / 2;
        }
        surface_area.push_back(area);
    }
}

float calc_area(Vec p1, Vec p2, Vec p3) {
    return (.5) * ((p1 - p3) / (p2 - p3)).mag();
}


void init() {    

    // read in ply file
    read_test((char*)FILE_NAME);

    internal_diagonal = sqrtf(powf(min_corner.x - max_corner.x, 2) + powf(min_corner.y - max_corner.y, 2) + powf(min_corner.z - max_corner.z, 2));

    calc_normals();
    calc_face_area();
    cam = { Vec((min_corner.x + max_corner.x) / 2 - internal_diagonal,(min_corner.y + max_corner.y) / 2,(min_corner.z + max_corner.z) / 2), Vec(1,0,0), Vec(0,1,0), internal_diagonal, internal_diagonal * (((float)HEIGHT) / WIDTH), internal_diagonal };
    li = { Vec(255, 255, 255), Vec(255, 0, 0), .25, .25, .5 };
}

float clamp(float i, float min, float max) {
    if (i < min) i = min;
    if (i > max) i = max;
    return i;
}

void display() {
    glClearColor(0, 0, 0, 1);
    glClear(GL_COLOR_BUFFER_BIT);

    for (int y = 0; y < HEIGHT; ++y)
    {
        cout << "row " << y << "/" << HEIGHT << endl;
        for (int x = 0; x < WIDTH; ++x)
        {
     
            Vec p = (cam.eye) + 
                    (cam.vd * cam.dist) + 
                    (cam.up * (y - HEIGHT / 2) * (cam.h / (HEIGHT / 2))) + 
                    ((cam.up / cam.vd) * (x - WIDTH / 2) * (cam.w / (WIDTH / 2)));
            
            Vec v = (p - cam.eye) / (p - cam.eye).mag();
            int closestFace = 0;
            Vec closestPoint;
            float closest_dist = INFINITY;
            for (int i = 0; i < num_elems; i++) {
                float t = ((Vec(*vlist[flist[i]->verts[0]]) - cam.eye) * surface_normals[i]) / (v * surface_normals[i]);
                Vec s = cam.eye + v * t;
                float area = calc_area(Vec(*vlist[flist[i]->verts[flist[i]->nverts - 1]]), Vec(*vlist[flist[i]->verts[0]]), s);
                for (int j = 0; j < flist[i]->nverts - 1; j++) {
                    area += calc_area(Vec(*vlist[flist[i]->verts[j]]), Vec(*vlist[flist[i]->verts[j + 1]]), s);
                }
                if (t > 0 && t < cam.eye.mag() + internal_diagonal && area <= surface_area[i] * (1 + fuzz)) {
                    if (t < closest_dist) {
                        closestFace = i;
                        closestPoint = s;
                        closest_dist = t;
                    }
                }
            }
            //std::cout << x << ", " << y << ": " << closest_dist << std::endl;
            if (closest_dist < INFINITY) {

                Vec L = cam.eye - closestPoint;
                L.normalize();
                Vec N = surface_normals[closestFace];
                Vec R = N * (2 * (N * L)) - L;

                Vec I = li.i_i * (li.k_d * (L * N) + li.k_s * powf(R * L, 10)) + li.i_a * li.k_a;
                
                pixels[y][x][0] = clamp(I.x, 0, 255) * 256 * 256 * 256;
                pixels[y][x][1] = clamp(I.y, 0, 255) * 256 * 256 * 256;
                pixels[y][x][2] = clamp(I.z, 0, 255) * 256 * 256 * 256;
            }
            else {
                pixels[y][x][0] = (0) * 256 * 256 * 256;
                pixels[y][x][1] = (0) * 256 * 256 * 256;
                pixels[y][x][2] = (0) * 256 * 256 * 256;

          
            }
        }
    }
    
    glDrawPixels(WIDTH, HEIGHT, GL_RGB, GL_UNSIGNED_INT, pixels);

	glFlush();

    
}

Vec rotateAbout(Vec v, Vec axis, float theta) {
    return axis * (axis * v) + ((axis / v) / axis) * cosf(theta) + (axis / v) * sinf(theta);
}

void KeyboardFunc(unsigned char Key, int x, int y)
{
    float step = .1;
    Vec axis = cam.up / cam.vd;
    axis.normalize();

    switch (Key)
    {
    case 'w':
        cout << "zooming in" << endl;
        cout << "moved from: " << cam.eye << endl;
        cam.eye = cam.eye + (cam.vd * step);
        cout << "move to: " << cam.eye << endl;
        break;
    case 'a':
        cout << "panning left" << endl;
        cout << "changed looking at from: " << cam.vd << endl;
        cam.vd = rotateAbout(cam.vd, cam.up, -0.1);
        cout << "changed looking at to: " << cam.vd << endl;
        break;
    case 's':
        cout << "zooming out" << endl;
        cout << "moved from: " << cam.eye << endl;
        cam.eye = cam.eye - (cam.vd * step);
        cout << "move to: " << cam.eye << endl;
        break;
    case 'd':
        cout << "panning right" << endl;
        cout << "changed looking at from: " << cam.vd << endl;
        cam.vd = rotateAbout(cam.vd, cam.up, 0.1);
        cout << "changed looking at to: " << cam.vd << endl;
        break;
    case 'q':
        cout << "panning up" << endl;
        cout << "changed looking at from: " << cam.vd << endl;
        cam.vd = rotateAbout(cam.vd, axis, -0.1);
        cout << "changed looking at to: " << cam.vd << endl;
        cam.up = rotateAbout(cam.up, axis, -0.1);
        break;
    case 'e':
        cout << "panning down" << endl;
        cout << "changed looking at from: " << cam.vd << endl;
        cam.vd = rotateAbout(cam.vd, axis, 0.1);
        cout << "changed looking at to: " << cam.vd << endl;
        cam.up = rotateAbout(cam.up, axis, 0.1);
        break;
    case ' ':
        cout << "saving to ./output.ppm" << endl;
        save_ppm("output.ppm");
        cout << "saved" << endl;
        break;
    default:
        break;
    };

    cam.vd.normalize();
    cam.up.normalize();
    glutPostRedisplay();
}

int main(int argc, char** argv)
{
	glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_SINGLE | GLUT_RGB);
	glutInitWindowPosition(50, 50);
	glutInitWindowSize(WIDTH, HEIGHT);
	glutCreateWindow("Assignment 4");
	init();
	glutDisplayFunc(display);
    glutKeyboardFunc(KeyboardFunc);
	glutMainLoop();

}

// Run program: Ctrl + F5 or Debug > Start Without Debugging menu
// Debug program: F5 or Debug > Start Debugging menu

// Tips for Getting Started: 
//   1. Use the Solution Explorer window to add/manage files
//   2. Use the Team Explorer window to connect to source control
//   3. Use the Output window to see build output and other messages
//   4. Use the Error List window to view errors
//   5. Go to Project > Add New Item to create new code files, or Project > Add Existing Item to add existing code files to the project
//   6. In the future, to open this project again, go to File > Open > Project and select the .sln file
