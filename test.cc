#include "Pythia8/Pythia.h"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/AreaDefinition.hh"
#include "fastjet/GhostedAreaSpec.hh"
#include "TH1F.h"
#include "TH2F.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TLegend.h"
#include <fstream>
#include <iostream>
#include <cmath>
#include <vector>
// salam gavin , matteo caccrai
// Define a class to represent a particle or constituent
// min jet pt 5/10
class Particle {
public:
    double pt;
    double eta;
    double phi;

    Particle(double _pt, double _eta, double _phi) : pt(_pt), eta(_eta), phi(_phi) {}
};

void saveHistogram2DAsROOT(TH2F &histogram, const char* filename) {
    // Create a canvas
    TCanvas canvas("canvas", "Histogram Canvas");

    // Draw the histogram on the canvas
    histogram.Draw("colz");  // "colz" option for a 2D color plot

    // Save the canvas and histogram in a ROOT file
    std::string fullFilename = std::string("Histograms") + "/" + std::string(filename);
    TFile outputFile(fullFilename.c_str(), "RECREATE");
    canvas.Write();
    histogram.Write();
    outputFile.Close();
}




void plotOverlayHistograms(TH1F& hist1, TH1F& hist2, const char* label1, const char* label2, const char* filename) {
    // Create a canvas
    TCanvas canvas("canvas", filename, 800, 600);

    // Draw the first histogram with blue color
    hist1.SetLineColor(kBlue);
    hist1.Draw();

    // Draw the second histogram with red color, overlaid on top of the first histogram
    hist2.SetLineColor(kRed);
    hist2.Draw("same");

    // Add a legend
    TLegend legend(0.88, 0.6, 0.98, 0.7);
    legend.AddEntry(&hist1, label1, "l");
    legend.AddEntry(&hist2, label2, "l");
    legend.Draw();

    // Save the canvas as a ROOT file
    std::string fullFilename = std::string("Histograms") + "/" + std::string(filename);
    TFile outputFile(fullFilename.c_str(), "RECREATE");
    canvas.Write();
    hist1.Write();
    hist2.Write();
    outputFile.Close();
}


// Function to calculate Euclidean distance between two points in the eta-phi plane
double euclideanDistance(double eta1, double phi1, double eta2, double phi2) {
    double deta = eta1 - eta2;
    double dphi = phi1 - phi2;
    while (dphi > M_PI) dphi -= 2 * M_PI;
    while (dphi <= -M_PI) dphi += 2 * M_PI;
    return sqrt(deta * deta + dphi * dphi);
}

int main() {
    // Create Pythia instance
    Pythia8::Pythia pythia_g;
    Pythia8::Pythia pythia_q;
    int event_size = 1000000;

    // Set up Pythia configuration
    pythia_g.readString("Beams:eCM = 13000."); // Center-of-mass energy
    pythia_g.readString("PhaseSpace:pTHatMin = 20.");
    pythia_g.readString("Tune:pp = 14");
    pythia_g.readString("PartonLevel:ISR = on");  // Enable initial-state radiation

    pythia_q.readString("Beams:eCM = 13000."); // Center-of-mass energy
    pythia_q.readString("PhaseSpace:pTHatMin = 20.");
    pythia_q.readString("Tune:pp = 14");
    pythia_q.readString("PartonLevel:ISR = on");  // Enable initial-state radiation

    pythia_g.readString("HardQCD:qqbar2gg = on");  // Turn on hard QCD processes
    pythia_g.readString("HardQCD:gg2gg = on");
    pythia_g.readString("HardQCD:gg2ggg = on");

    pythia_q.readString("HardQCD:gg2qqbar = on");
    pythia_q.readString("HardQCD:qq2qq = on");
    

    // Initialize Pythia
    if (!pythia_g.init()) {
        std::cerr << "Initialization failed!" << std::endl;
        return 1;
    }

    // Define jet radius
    double jet_R = 0.5;
    fastjet::Strategy               strategy = fastjet::Best;
    fastjet::RecombinationScheme    recombScheme = fastjet::E_scheme;

    // Initialize FastJet jet definition
    fastjet::JetDefinition jet_def(fastjet::antikt_algorithm, jet_R,
                                      recombScheme, strategy);

    // Initialize FastJet area definition
    fastjet::AreaDefinition area_def(fastjet::active_area_explicit_ghosts, fastjet::GhostedAreaSpec(5.0));
    
    // TCanvas canvas("canvas", "Number of Jets", 800, 600);
    // Open a file to write CSV data
    std::ofstream outFile("event_properties_2M.csv");

    // Write header line with column names
    outFile << "no_of_jets,lead_jet_pt,lead_jet_eta,lead_jet_phi,lead_jet_mass,lead_jet_energy,lead_jet_multiplicity,lead_jet_m/pt,avg_pt,var_pt,avg_distance,var_distance,var_phi,jet_width_1,jet_width_2,sublead_jet_pt,sublead_jet_eta,sublead_jet_phi,target" << std::endl;
    
    TH2F gluon_phi_eta_hist("hist_PhiVsEta", "Phi vs Eta Distribution", 100, -5, 5, 100, -3.14, 3.14);
    TH1F gluon_no_of_jets("hist_gluon_no_of_jets", "Number of Jets per Event", 50, 0, 50);
    TH1F gluon_lead_jet_pt("hist_gluon_lead_jet_pt", "Leading Jet pT", 200, 0, 100); //change
    TH1F gluon_lead_jet_energy("hist_gluon_lead_jet_energy", "Leading Jet Energy", 200, 0, 500);
    TH1F gluon_lead_jet_multiplicity("hist_gluon_lead_jet_multiplicity", "Leading Jet Multiplicity", 100, 0, 100); //change
    TH1F gluon_lead_jet_mass_over_pt("hist_gluon_lead_jet_mass_over_pt", "Leading Jet Mass/pT", 200, 0, 2); //change
    TH1F gluon_avg_pt("hist_gluon_avg_pt", "Average pT of Jets", 200, 0, 20);
    TH1F gluon_var_pt("hist_gluon_var_pt", "Variance of pT of Jets", 200, 0, 100);
    TH1F gluon_avg_distance("hist_gluon_avg_distance", "Average Distance of Jets", 200, 0, 1);
    TH1F gluon_jet_width_1("hist_gluon_jet_width_1", "Jet Width 1 of Jets", 200, 0, 1);  // chnage
    TH1F gluon_sublead_jet_pt("hist_gluon_sublead_jet_pt", "Subleading Jet pT", 200, 0, 100); //change

    
    // Loop over events
    for (int iEvent = 0; iEvent < event_size; ++iEvent) {
        // Generate event
        if (!pythia_g.next()) continue;

        // Convert Pythia particles to FastJet PseudoJets
        std::vector<fastjet::PseudoJet> particles;
        for (int i = 0; i < pythia_g.event.size(); ++i) {
            if (pythia_g.event[i].isFinal() && pythia_g.event[i].isVisible()) {
                particles.push_back(fastjet::PseudoJet(pythia_g.event[i].px(), pythia_g.event[i].py(), pythia_g.event[i].pz(), pythia_g.event[i].e()));
                gluon_phi_eta_hist.Fill(pythia_g.event[i].eta(), pythia_g.event[i].phi());
            }
        }

        // Cluster particles into jets using FastJet
        fastjet::ClusterSequence cs(particles, jet_def);
        std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets(5));
        // Process each jet

        if (!jets.empty()) {
            double lead_jet_pt = jets[0].pt();
            double lead_jet_eta = jets[0].eta();
            double lead_jet_phi = jets[0].phi();
            double lead_jet_mass = jets[0].m();
            double lead_jet_energy = jets[0].e();
            double lead_jet_multiplicity = jets[0].constituents().size();
            double lead_jet_mass_over_pt = lead_jet_mass / lead_jet_pt;

            // Calculate average pt and variance
            double avg_pt = 0.0;
            double var_pt = 0.0;
            double avg_distance = 0.0;
            double var_distance = 0.0;
            double sum_weighted_pt_1 = 0.0;
            double sum_weighted_pt_2 = 0.0;
            double sum_phi_sq = 0.0;
            double sum_phi = 0.0;

            for (const auto& constituent : jets[0].constituents()) {
                double pt = constituent.pt();
                avg_pt += pt;
                var_pt += pt * pt;
                double distance = euclideanDistance(lead_jet_eta, lead_jet_phi, constituent.eta(), constituent.phi());
                avg_distance += distance;
                var_distance += distance * distance;
                sum_weighted_pt_1 += pt * distance;
                sum_weighted_pt_2 += pt * std::sqrt(distance);
                sum_phi += constituent.phi();
                sum_phi_sq += pow(constituent.phi(),2); 
            }
            avg_pt /= lead_jet_multiplicity;
            var_pt = var_pt / lead_jet_multiplicity - avg_pt * avg_pt;
            avg_distance /= lead_jet_multiplicity;
            var_distance = var_distance / lead_jet_multiplicity - avg_distance * avg_distance;
            double jet_width_1 = sum_weighted_pt_1 / (avg_pt*lead_jet_multiplicity);
            double jet_width_2 = sum_weighted_pt_2 / (avg_pt*lead_jet_multiplicity);
            double var_phi = sum_phi_sq / lead_jet_multiplicity - pow(sum_phi / lead_jet_multiplicity, 2);
            
            // Calculate properties for subleading jet
            double sublead_jet_pt = 0.0;
            double sublead_jet_eta = 0.0;
            double sublead_jet_phi = 0.0;
            if (jets.size() > 1) {
                sublead_jet_pt = jets[1].pt();
                sublead_jet_eta = jets[1].eta();
                sublead_jet_phi = jets[1].phi();
            }

            // Write properties to CSV file
            outFile << jets.size() << ","
                    << lead_jet_pt << ","
                    << lead_jet_eta << ","
                    << lead_jet_phi << ","
                    << lead_jet_mass << ","
                    << lead_jet_energy << ","
                    << lead_jet_multiplicity << ","
                    << lead_jet_mass_over_pt << ","
                    << avg_pt << ","
                    << var_pt << ","
                    << avg_distance << ","
                    << var_distance << ","
                    << var_phi << ","
                    << jet_width_1 << ","
                    << jet_width_2 << ","
                    << sublead_jet_pt << ","
                    << sublead_jet_eta << ","
                    << sublead_jet_phi << ","
                    << 1 << std::endl;

                    gluon_no_of_jets.Fill(jets.size());
                    gluon_lead_jet_pt.Fill(lead_jet_pt);
                    gluon_lead_jet_energy.Fill(lead_jet_energy);
                    gluon_lead_jet_multiplicity.Fill(lead_jet_multiplicity);
                    gluon_lead_jet_mass_over_pt.Fill(lead_jet_mass_over_pt);
                    gluon_avg_pt.Fill(avg_pt);
                    gluon_var_pt.Fill(var_pt);
                    gluon_avg_distance.Fill(avg_distance);
                    gluon_jet_width_1.Fill(jet_width_1);
                    gluon_sublead_jet_pt.Fill(sublead_jet_pt);

        } else {
            // If no jets found, write placeholder values to CSV file
            outFile << "0,N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A,1" << std::endl;
        }
        
        
    }



    std::ofstream statsFile("pythia_g_stats.txt");

    // Redirect the output of pythia.stat() to the file stream
    std::streambuf* orig_cout = std::cout.rdbuf();  // Save original cout buffer
    std::cout.rdbuf(statsFile.rdbuf());             // Redirect cout to the file

    // Call pythia.stat() to print statistics
    pythia_g.stat();

    // Restore the original cout buffer
    std::cout.rdbuf(orig_cout);

    // Close the file
    statsFile.close();

    // saveHistogramAsPNG(gluon_pT, "gluon_pT");
    // saveHistogramAsPNG(gluon_jets, "gluon_jets");
    // saveHistogramAsPNG(gluon_jet_multiplicity, "gluon_jet_multiplicity");



    if (!pythia_q.init()) {
        std::cerr << "Initialization failed!" << std::endl;
        return 1;
    }

    TH2F quark_phi_eta_hist("hist_PhiVsEta", "Phi vs Eta Distribution", 100, -5, 5, 100, -3.14, 3.14);
    TH1F quark_no_of_jets("hist_quark_no_of_jets", "Number of Jets per Event", 50, 0, 50);
    TH1F quark_lead_jet_pt("hist_quark_lead_jet_pt", "Leading Jet pT", 200, 0, 100); //change
    TH1F quark_lead_jet_energy("hist_quark_lead_jet_energy", "Leading Jet Energy", 200, 0, 500);
    TH1F quark_lead_jet_multiplicity("hist_quark_lead_jet_multiplicity", "Leading Jet Multiplicity", 100, 0, 100); //change
    TH1F quark_lead_jet_mass_over_pt("hist_quark_lead_jet_mass_over_pt", "Leading Jet Mass/pT", 200, 0, 2); //change
    TH1F quark_avg_pt("hist_quark_avg_pt", "Average pT of Jets", 200, 0, 20);
    TH1F quark_var_pt("hist_quark_var_pt", "Variance of pT of Jets", 200, 0, 100);
    TH1F quark_avg_distance("hist_quark_avg_distance", "Average Distance of Jets", 200, 0, 1);
    TH1F quark_jet_width_1("hist_quark_jet_width_1", "Jet Width 1 of Jets", 200, 0, 1);  // change
    TH1F quark_sublead_jet_pt("hist_quark_sublead_jet_pt", "Subleading Jet pT", 200, 0, 100); // change


    // Loop over events
    for (int iEvent = 0; iEvent < event_size; ++iEvent) {
        // Generate event
        if (!pythia_q.next()) continue;

        // Convert Pythia particles to FastJet PseudoJets
        std::vector<fastjet::PseudoJet> particles;
        for (int i = 0; i < pythia_q.event.size(); ++i) {
            if (pythia_q.event[i].isFinal() && pythia_q.event[i].isVisible()) {
                particles.push_back(fastjet::PseudoJet(pythia_q.event[i].px(), pythia_q.event[i].py(), pythia_q.event[i].pz(), pythia_q.event[i].e()));
                quark_phi_eta_hist.Fill(pythia_q.event[i].eta(), pythia_q.event[i].phi());
            }
        }

        // Cluster particles into jets using FastJet
        fastjet::ClusterSequence cs(particles, jet_def);
        std::vector<fastjet::PseudoJet> jets = fastjet::sorted_by_pt(cs.inclusive_jets(5));
        // Process each jet
    

        if (!jets.empty()) {
            double lead_jet_pt = jets[0].pt();
            double lead_jet_eta = jets[0].eta();
            double lead_jet_phi = jets[0].phi();
            double lead_jet_mass = jets[0].m();
            double lead_jet_energy = jets[0].e();
            double lead_jet_multiplicity = jets[0].constituents().size();
            double lead_jet_mass_over_pt = lead_jet_mass / lead_jet_pt;

            // Calculate average pt and variance
            double avg_pt = 0.0;
            double var_pt = 0.0;
            double avg_distance = 0.0;
            double var_distance = 0.0;
            double sum_weighted_pt_1 = 0.0;
            double sum_weighted_pt_2 = 0.0;
            double sum_phi_sq = 0.0;
            double sum_phi = 0.0;

            for (const auto& constituent : jets[0].constituents()) {
                double pt = constituent.pt();
                avg_pt += pt;
                var_pt += pt * pt;
                double distance = euclideanDistance(lead_jet_eta, lead_jet_phi, constituent.eta(), constituent.phi());
                avg_distance += distance;
                var_distance += distance * distance;
                sum_weighted_pt_1 += pt * distance;
                sum_weighted_pt_2 += pt * std::sqrt(distance);
                sum_phi += constituent.phi();
                sum_phi_sq += pow(constituent.phi(),2); 
            }
            avg_pt /= lead_jet_multiplicity;
            var_pt = var_pt / lead_jet_multiplicity - avg_pt * avg_pt;
            avg_distance /= lead_jet_multiplicity;
            var_distance = var_distance / lead_jet_multiplicity - avg_distance * avg_distance;
            double jet_width_1 = sum_weighted_pt_1 / (avg_pt*lead_jet_multiplicity);
            double jet_width_2 = sum_weighted_pt_2 / (avg_pt*lead_jet_multiplicity);
            double var_phi = sum_phi_sq / lead_jet_multiplicity - pow(sum_phi / lead_jet_multiplicity, 2);
            
            // Calculate properties for subleading jet
            double sublead_jet_pt = 0.0;
            double sublead_jet_eta = 0.0;
            double sublead_jet_phi = 0.0;
            if (jets.size() > 1) {
                sublead_jet_pt = jets[1].pt();
                sublead_jet_eta = jets[1].eta();
                sublead_jet_phi = jets[1].phi();
            }

            // Write properties to CSV file
            outFile << jets.size() << ","
                    << lead_jet_pt << ","
                    << lead_jet_eta << ","
                    << lead_jet_phi << ","
                    << lead_jet_mass << ","
                    << lead_jet_energy << ","
                    << lead_jet_multiplicity << ","
                    << lead_jet_mass_over_pt << ","
                    << avg_pt << ","
                    << var_pt << ","
                    << avg_distance << ","
                    << var_distance << ","
                    << var_phi << ","
                    << jet_width_1 << ","
                    << jet_width_2 << ","
                    << sublead_jet_pt << ","
                    << sublead_jet_eta << ","
                    << sublead_jet_phi << ","
                    << 0 << std::endl;

                    quark_no_of_jets.Fill(jets.size());
                    quark_lead_jet_pt.Fill(lead_jet_pt);
                    quark_lead_jet_energy.Fill(lead_jet_energy);
                    quark_lead_jet_multiplicity.Fill(lead_jet_multiplicity);
                    quark_lead_jet_mass_over_pt.Fill(lead_jet_mass_over_pt);
                    quark_avg_pt.Fill(avg_pt);
                    quark_var_pt.Fill(var_pt);
                    quark_avg_distance.Fill(avg_distance);
                    quark_jet_width_1.Fill(jet_width_1);
                    quark_sublead_jet_pt.Fill(sublead_jet_pt);


        } else {
            // If no jets found, write placeholder values to CSV file
            outFile << "0,N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A,N/A,0" << std::endl;
        }

        
    
    }

    std::ofstream statsFile_q("pythia_q_stats.txt");

    // Redirect the output of pythia.stat() to the file stream
    std::streambuf* orig_cout_q = std::cout.rdbuf();  // Save original cout buffer
    std::cout.rdbuf(statsFile_q.rdbuf());             // Redirect cout to the file

    // Call pythia.stat() to print statistics
    pythia_q.stat();

    // Restore the original cout buffer
    std::cout.rdbuf(orig_cout_q);

    // Close the file
    statsFile.close();

    outFile.close();

        // Save histograms as PNG files
    // saveHistogramAsPNG(quark_pT, "quark_pT");
    // saveHistogramAsPNG(quark_jets, "quark_jets");
    // saveHistogramAsPNG(quark_jet_multiplicity, "quark_jet_multiplicity");

    // Save histograms
    plotOverlayHistograms(gluon_no_of_jets, quark_no_of_jets, "Gluon", "Quark", "no_of_jets_overlay.root");
    plotOverlayHistograms(gluon_lead_jet_pt, quark_lead_jet_pt, "Gluon", "Quark", "lead_jet_pt_overlay.root"); //
    plotOverlayHistograms(gluon_lead_jet_energy, quark_lead_jet_energy, "Gluon", "Quark", "lead_jet_energy_overlay.root");
    plotOverlayHistograms(gluon_lead_jet_multiplicity, quark_lead_jet_multiplicity, "Gluon", "Quark", "lead_jet_multiplicity_overlay.root"); //
    plotOverlayHistograms(gluon_lead_jet_mass_over_pt, quark_lead_jet_mass_over_pt, "Gluon", "Quark", "lead_jet_mass_over_pt_overlay.root"); //
    plotOverlayHistograms(gluon_avg_pt, quark_avg_pt, "Gluon", "Quark", "avg_pt_overlay.root");
    plotOverlayHistograms(gluon_var_pt, quark_var_pt, "Gluon", "Quark", "var_pt_overlay.root");
    plotOverlayHistograms(gluon_avg_distance, quark_avg_distance, "Gluon", "Quark", "avg_distance_overlay.root");
    plotOverlayHistograms(gluon_jet_width_1, quark_jet_width_1, "Gluon", "Quark", "jet_width_1_overlay.root"); //
    plotOverlayHistograms(gluon_sublead_jet_pt, quark_sublead_jet_pt, "Gluon", "Quark", "sublead_jet_pt_overlay.root"); //

    saveHistogram2DAsROOT(quark_phi_eta_hist, "quark_phi_eta_hist.root");
    saveHistogram2DAsROOT(gluon_phi_eta_hist, "gluon_phi_eta_hist.root");


    // Save Pythia statistics to files
    

    


    return 0;
}