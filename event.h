#include "TLorentzVector.h"
#include <deque>
#include <fstream>

std::vector<std::vector<double>> LUND_OUTPUT;

struct Particle
{
    TLorentzVector p;
    int id;
    double mass;
};

class Event
{
    private:
        double x, y, z, cross_section, beam_energy;
        std::deque<Particle> bunch;
        int beam_polarization;
        bool cm_system;
    public:
        TVector3 get_coordinates()
        {
            TVector3 beta(x, y, z);
            return beta;
        }

        Event();
        Event(TVector3&);
        Event(TVector3&, double&);

        Event(const Event& buff) : x(buff.x), y(buff.y), z(buff.z), cross_section(buff.cross_section), beam_energy(buff.beam_energy),
        bunch(buff.bunch), beam_polarization(buff.beam_polarization), cm_system(buff.cm_system) {}

        double get_section(){return cross_section;}
        void set_beam(double&, int&);
        void set_coordinates(double&, double&);
        void set_section(double&);
        void add_particle(Particle&);
        void clear_event(){bunch.clear();}
        void print_vector();
        void print_lund(std::string&);
        void set_cm_system();
        void set_lab_system();
        void cm_to_lab(double&, double&);
        void lab_to_cm(double&, double&);
        void Z_rotate_random();


};

Event::Event()
{
    x = 0; y = 0; z = 0;
    beam_energy = 6.565;
    cross_section = 0.001;
    beam_polarization = 0;
    cm_system = true;
}

Event::Event(TVector3& beta)
{
    x = beta.X(); y = beta.Y(); z = beta.Z();
    beam_energy = 6.565;
    cross_section = 0.001;
    beam_polarization = 0;
    cm_system = true;
}

Event::Event(TVector3& beta, double& E)
{
    x = beta.X(); y = beta.Y(); z = beta.Z();
    beam_energy = E;
    cross_section = 0.001;
    beam_polarization = 0;
    cm_system = true;
}

void Event::set_beam(double& E, int& h)
{
    beam_energy = E;
    beam_polarization = h;
}

void Event::set_coordinates(double& R, double& L)
{
    x = 2*R;
    y = 2*R;

    while(x*x + y*y > R*R)
    {
        x = fRand(-R, R);
        y = fRand(-R, R);
    }
    z = fRand(-L/2-2.5, L/2-2.5);
}

void Event::set_section(double& S)
{
    cross_section = S;
}

void Event::add_particle(Particle& buff)
{
    bunch.push_back(buff);
}

void Event::print_vector()
{
    std::vector<double> buff;
    buff.push_back(bunch.size());
    buff.push_back(1); buff.push_back(1); buff.push_back(0);
    buff.push_back(beam_polarization); buff.push_back(11);
    buff.push_back(beam_energy); buff.push_back(2212);
    buff.push_back(0); buff.push_back(cross_section);

    LUND_OUTPUT.push_back(buff); buff.clear();

    for(long unsigned int i = 0; i < bunch.size(); i++)
    {
        buff.push_back(i+1); buff.push_back(0); buff.push_back(1);
        buff.push_back(bunch[i].id); buff.push_back(0); buff.push_back(0);
        buff.push_back((bunch[i].p).Px()); buff.push_back((bunch[i].p).Py());
        buff.push_back((bunch[i].p).Pz()); buff.push_back((bunch[i].p).E());
        buff.push_back(bunch[i].mass);
        buff.push_back(x); buff.push_back(y); buff.push_back(z);

        LUND_OUTPUT.push_back(buff); buff.clear();
    }
}

void Event::print_lund(std::string& Path)
{
    std::ofstream File;

    File.open(Path,std::fstream::in | std::fstream::out | std::fstream::app);

    File << bunch.size() << "\t1\t1\t0\t" << beam_polarization << "\t11\t" << beam_energy << "\t2212\t0\t" << cross_section << std::endl;

    for(long unsigned int i = 0; i < bunch.size(); i++)
    {
        File << i+1 << "\t0\t1\t" << bunch[i].id << "\t0\t0\t" << (bunch[i].p).Px() << "\t" <<  (bunch[i].p).Py() << "\t"
        << (bunch[i].p).Pz() << "\t" << (bunch[i].p).E() << "\t" << bunch[i].mass << "\t" << x << "\t" << y << "\t" << z << std::endl;
    }

    File.close();
}

void Event::set_cm_system()
{
    cm_system = true;
}

void Event::set_lab_system()
{
    cm_system = false;
}

void Event::cm_to_lab(double& W, double& Q2)
{
    if(cm_system)
    {
        TVector3 beta; double nu;
        nu =  (W*W + Q2 - m_p*m_p)/(2*m_p);
        beta.SetXYZ(0., 0., sqrt(nu*nu + Q2)/(nu + m_p));
        double ang = fRand(0, 2*M_PI);

        for(long unsigned int i = 0; i < bunch.size(); i++)
        {
            (bunch[i].p).Boost(beta);
            (bunch[i].p).RotateY(-acos((Q2 + 2*beam_energy*nu)/(2*beam_energy*sqrt(nu*nu + Q2))));
            (bunch[i].p).RotateZ(ang);
        }
        cm_system = false;
    }
}

void Event::Z_rotate_random()
{
    double ang = fRand(0, 2*M_PI);

    for(long unsigned int i = 0; i < bunch.size(); i++)
    {
        (bunch[i].p).RotateZ(ang);
    }
}

void Event::lab_to_cm(double& W, double& Q2)
{
    if(cm_system == false)
    {
        TVector3 beta; double nu, ang1, ang2;
        TLorentzVector e_i, q; e_i.SetPxPyPzE(0, 0, beam_energy, beam_energy);

        for(long unsigned int i = 0; i < bunch.size(); i++)
        {
            if(bunch[i].id == 11)
            {
                ang1 = (bunch[i].p).Phi();
                q = e_i - bunch[i].p;
                ang2 = q.Theta();
            }
        }

        nu =  (W*W + Q2 - m_p*m_p)/(2*m_p);
        beta.SetXYZ(0., 0., -sqrt(nu*nu + Q2)/(nu + m_p));

        for(long unsigned int i = 0; i < bunch.size(); i++)
        {
            (bunch[i].p).RotateZ(-ang1);
            (bunch[i].p).RotateY(ang2);
            (bunch[i].p).Boost(beta);
        }
        cm_system = true;
    }
}
