#include "Pythia8/Pythia.h"
#include "Pythia8Plugins/HepMC2.h"
#include "HepMC/GenEvent.h"
#include "HepMC/IO_GenEvent.h"
#include "TTree.h"
#include "TFile.h"
#include "TLorentzVector.h"
#include <iostream>
#include <vector>
#include <cmath>
#include "TH1F.h"

using namespace Pythia8;
Int_t pdgid, nTracks; // nova variável para o PDG ID

int main(int argc, char *argv[])
{
  // Configurar e inicializar Pythia
  std::string arg1, arg2, arg3, arg4, nEvents;
  Pythia pythia;

  if (argc > 1)
  {
    arg1 = argv[1];
    if (arg1 == "MPICheck_off")
    {
      pythia.readString("Diffraction:sampleType = 3");
      std::cout << "Diffraction:sampleType = 3" << std::endl;
    }
    else if (arg1 == "MPICheck_on")
    {
      pythia.readString("Diffraction:sampleType = 4");
      std::cout << "Diffraction:sampleType = 4" << std::endl;
    }
    else
    {
      std::cout << "Usage: ./diff_qqbar MPICheck_on/off bbbar/ttbar/ccbar pomset5/pomset12/pomset13 nEvents" << std::endl;
      return 1;
    }
  }

  if (argc > 2)
  {
    arg2 = argv[2];
    if (arg2 == "bbbar")
    {
      pythia.readString("HardQCD:hardbbbar = on");
      std::cout << "HardQCD:qqbar2bbbar = on, qq + gg" << std::endl;
    }
    else if (arg2 == "ttbar")
    {
      pythia.readString("Top:qqbar2ttbar = on");
      std::cout << "Top:qqbar2ttbar = on, qq + gg" << std::endl;
    }
    else if (arg2 == "ccbar")
    {
      pythia.readString("HardQCD:hardccbar = on");
      std::cout << "HardQCD:qqbar2ccbar = on, qq + gg" << std::endl;
    }
    else
    {
      std::cout << "Usage: ./diff_qqbar MPICheck_on/off bbbar/ttbar/ccbar pomset5/pomset12/pomset13 nEvents" << std::endl;
      return 1;
    }
  }

  if (argc > 3)
  {
    arg3 = argv[3];
    if (arg3 == "pomset5")
    {
      pythia.readString("PDF:PomSet = 5");
      pythia.readString("SigmaDiffractive:PomFlux = 7");
      std::cout << "PDF:PomSet = 5 - h1jets2007 - PomFlux H1FitA" << std::endl;
    }
    else if (arg3 == "pomset12")
    {
      pythia.readString("PDF:PomSet = 12");
      pythia.readString("SigmaDiffractive:PomFlux = 6");
      std::cout << "PDF:PomSet = 12 - GKG18-DPDF-FitA - PomFlux H1FitA" << std::endl;
    }
    else if (arg3 == "pomset13")
    {
      pythia.readString("PDF:PomSet = 13");
      pythia.readString("SigmaDiffractive:PomFlux = 7");
      std::cout << "PDF:PomSet = 13 -GKG18-DPDF-FitB- PomFlux H1Fit B" << std::endl;
    }
    else
    {
      std::cout << "Usage: ./diff_qqbar MPICheck_on/off bbbar/ttbar/ccbar pomset5/pomset12/pomset13 nEvents" << std::endl;
      return 1;
    }
  }
  if (argc > 4)
  {
    std::string arg4 = argv[4]; // Converte arg4 para std::string
    if (arg4 == "long")
    {
      pythia.readString("Main:numberOfEvents = 100000");
      std::cout << "Main:numberOfEvents = 100000" << std::endl;
    }
    else if (arg4 == "short")
    {
      pythia.readString("Main:numberOfEvents = 1000");
      std::cout << "Main:numberOfEvents = 1000" << std::endl;
    }
    else if (arg4 == "medium")
    {
      pythia.readString("Main:numberOfEvents = 10000");
      std::cout << "Main:numberOfEvents = 10000" << std::endl;
    }
  }

  pythia.readFile("config_pp_qqbar_diffraction.cmnd");
  pythia.init();
  // Inicializar o escritor HepMC2
  HepMC::IO_GenEvent ascii_io("output.hepmc", std::ios::out);
  // Conversor de eventos Pythia8 para HepMC2
  HepMC::Pythia8ToHepMC toHepMC;
  nTracks = 0; // inicializar o contador de tracks
  // Inicializar o arquivo e a árvore do ROOT

  std::stringstream ss;

  ss << "output_ap_ttbar_tree_" << arg1 << ".root";

  TFile f(ss.str().c_str(), "RECREATE");

  TTree t("t", "State final tree");

  // Definir as variáveis que serão preenchidas na árvore
  Double_t px, py, pz, E, eta;
  Double_t b1px, b1py, b1pz, b1E, b1m;
  Double_t b2px, b2py, b2pz, b2E, b2m;
  Double_t mb1b2;
  Double_t protonpx, protonpy, protonpz, protonE, protonm, protoneta;

  TH1F *h_ntracks = new TH1F("ntracks", "ntracks", 500, 0, 500);
  TH1F *h_ntracks_sideA = new TH1F("ntracks_sideA", "ntracks_sideA", 500, 0, 500);
  TH1F *h_ntracks_sideB = new TH1F("ntracks_sideB", "ntracks_sideB", 500, 0, 500);
  TH1F Mb1b2("h_x", "", 100, -1e-4, 1e-4);

  double xi = 0;
  // Int_t pdgid, nTracks = 0; // nova variável para o PDG ID
  t.Branch("px", &px, "px/D");
  t.Branch("py", &py, "py/D");
  t.Branch("pz", &pz, "pz/D");
  t.Branch("E", &E, "E/D");
  t.Branch("eta", &eta, "eta/D");
  t.Branch("pdgid", &pdgid, "pdgid/I"); // nova branch para o PDG ID
  t.Branch("nTracks", &nTracks, "nTracks/I");

  t.Branch("xi", &xi, "xi/D");

  t.Branch("b1px", &b1px, "b1px/D");
  t.Branch("b1py", &b1py, "b1py/D");
  t.Branch("b1pz", &b1pz, "b1pz/D");
  t.Branch("b1E", &b1E, "b1E/D");
  t.Branch("b1m", &b1m, "b1m/D");

  t.Branch("b2px", &b2px, "b2px/D");
  t.Branch("b2py", &b2py, "b2py/D");
  t.Branch("b2pz", &b2pz, "b2pz/D");
  t.Branch("b2E", &b2E, "b2E/D");
  t.Branch("b2m", &b2m, "b2m/D");

  t.Branch("protonpx", &protonpx, "protonpx/D");
  t.Branch("protonpy", &protonpy, "protonpy/D");
  t.Branch("protonpz", &protonpz, "protonpz/D");
  t.Branch("protonE", &protonE, "protonE/D");
  t.Branch("protonm", &protonm, "protonm/D");
  t.Branch("protoneta", &protoneta, "protoneta/D");

  t.Branch("mb1b2", &mb1b2, "mb1b2/D");

  // Laço de eventos
  for (int iEvent = 0; iEvent < pythia.mode("Main:numberOfEvents"); ++iEvent)
  {
    // Gerar evento
    if (!pythia.next())
      continue;

    // Converter evento Pythia8 para HepMC2
    HepMC::GenEvent *hepmcEvent = new HepMC::GenEvent();
    toHepMC.fill_next_event(pythia, hepmcEvent);

    // Escrever evento HepMC2 no arquivo
    ascii_io << hepmcEvent;

    bool central = false;

    int nTracks_sideA = 0;
    int nTracks_sideB = 0;

    for (HepMC::GenEvent::particle_const_iterator p = hepmcEvent->particles_begin(); p != hepmcEvent->particles_end(); ++p)
    {
      if ((*p)->status() == 1)
      { // partícula no estado final
        px = (*p)->momentum().px();
        py = (*p)->momentum().py();
        pz = (*p)->momentum().pz();
        E = (*p)->momentum().e();
        eta = (*p)->momentum().eta();
        pdgid = (*p)->pdg_id();

        if (fabs(eta) < 2.5)
        {
          central = true;
        }

        if (pz < 2000 and !central)
        {
          nTracks_sideA++;
        }

        if (pz > -2000 and !central)
        {
          nTracks_sideB++;
        }
        nTracks++;
      }
    }
    h_ntracks->Fill(nTracks);
    h_ntracks_sideA->Fill(nTracks_sideA);
    h_ntracks_sideB->Fill(nTracks_sideB);

    for (HepMC::GenEvent::particle_const_iterator p = hepmcEvent->particles_begin(); p != hepmcEvent->particles_end(); ++p)
    {
      if ((*p)->pdg_id() == 5 && (*p)->status() == 23)
      { // quark b sem hadronizar
        b1px = (*p)->momentum().px();
        b1py = (*p)->momentum().py();
        b1pz = (*p)->momentum().pz();
        b1E = (*p)->momentum().e();
        b1m = (*p)->momentum().m();
      }
    }
    for (HepMC::GenEvent::particle_const_iterator p = hepmcEvent->particles_begin(); p != hepmcEvent->particles_end(); ++p)
    {
      if ((*p)->pdg_id() == -5 && (*p)->status() == 23)
      { // quark b sem hadronizar
        b2px = (*p)->momentum().px();
        b2py = (*p)->momentum().py();
        b2pz = (*p)->momentum().pz();
        b2E = (*p)->momentum().e();
        b2m = (*p)->momentum().m();
      }
    }

    TLorentzVector b1, b2, btot, protonp4;
    b1.SetPxPyPzE(b1px, b1py, b1pz, b1E);
    b2.SetPxPyPzE(b2px, b2py, b2pz, b2E);
    btot = b1 + b2;

    // Invariant mass of the b-quark pair
    mb1b2 = btot.M();
    
    // std::cout << "Invariant mass of the b-quark pair: " << btot.M() << std::endl;

    double beam = 6500;
    for (HepMC::GenEvent::particle_const_iterator p = hepmcEvent->particles_begin(); p != hepmcEvent->particles_end(); ++p)
    {
      if ((*p)->status() == 1 and (*p)->pdg_id() == 2212)
      {                                              // partícula no estado final, é um próton e tem momento na direção Z maior que 5100 GeV
        xi = 1 - fabs((*p)->momentum().pz()) / beam; // calcular a variável xi
        protonp4.SetPxPyPzE((*p)->momentum().px(), (*p)->momentum().py(), (*p)->momentum().pz(), (*p)->momentum().e());
        //    if (protonp4.Pz() < 0) {
        protonpx = protonp4.Px();
        protonpy = protonp4.Py();
        protonpz = protonp4.Pz();
        protonE = protonp4.E();
        protonm = protonp4.M();
        protoneta = protonp4.Eta();
        //  }
      }
    }

    delete hepmcEvent;
    t.Fill();
  }
  h_ntracks->Write();
  h_ntracks_sideA->Write();
  h_ntracks_sideB->Write();

  t.Write();
  f.Close();

  // Imprimir estatísticas finais
  pythia.stat();
  return 0;
}
