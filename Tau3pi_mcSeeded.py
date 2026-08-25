##################################################################################################
#####       Future Circular Collider Software - Tau-lepton Mass Peak Reconstruction          #####
#####                                  Steering file                                         #####
#####              Base file provided by Willy Weber and Aurelien Martens                    #####
##################################################################################################

import ROOT

processList = {
    "p8_ee_Ztautau_ecm91": {}
}

##################################################################################################
# Load custom C-functions
##################################################################################################
includePaths = ["analyzers_Tau3Pi.h"]

inputDir    = "/eos/experiment/fcc/ee/generation/DelphesEvents/winter2023/IDEA" # Whole sample
outputDir   = "result"

##################################################################################################
# Main class of the Analysis
##################################################################################################
class RDFanalysis():
    output_var = []
    
    def analysers(dframe):
        ##################################################################################################
        # Set up MC_history collections
        ##################################################################################################
        dframe0 = (
            dframe
            # UNCOMMENT TO USE A SMALLER SUBSET OF THE SAMPLE (FAST CODE TESTING). WHOLE SAMPLE: 100M EVENTS
            #.Filter("0 < rdfentry_ && rdfentry_<= 1000000")
            
            .Alias("Parents", "Particle#0.index") # parents
            .Alias("Daughters", "Particle#1.index") # daughters
            
            # MC indices of the decay
            # tau- (PDG = 15) -> pi- (PDG = -211) pi+ (PDG = 211) pi- (PDG = -211) nu_tau (16)
            # exclusive decay to remove events with additionnal pi0s, photons, other particles
            # Retrieves a vector of integers which correspond to indices in the
            # Particle block
            # vector[0] = the mother, and then the daughters in the order
            # specified, i.e. here [1] = the pi-, [2] = the pi+, [3] = the pi-,
            # [4] = the nu_tau
            #
            # Boolean arguments:
            #   1st: `stableDaughters`, when set to true, the daughters specified
            #        in the list are looked for among the final, stable
            #        particles that come out from the mother, i.e. the decay
            #        tree is explored recursively if needed.
            #   2nd: `chargeConjugateMother`
            #   3rd: `chargeConjugateDaughters`
            #   4th: `inclusiveDecay`, when set to false, if a mother is found,
            #        that decays into the particles specified in the list plus
            #        other particle(s), this decay is not selected.
            # If the event contains more than one such decays, only the first
            # one is kept.
            .Define("MCtau3pinu_indices",
                    "MCParticle::get_indices(15, {-211,211,-211,16}, true, false, false, false) (Particle, Daughters)")
            .Define("MCtau3pinu_indices_CC",
                    "MCParticle::get_indices(-15, {211,-211,211,-16}, true, false, false, false) (Particle, Daughters)")
            .Filter("MCtau3pinu_indices.size() > 0 && MCtau3pinu_indices_CC.size() > 0")

            .Define("MC_tau", "Particle.at(MCtau3pinu_indices[0])")
            .Define("MC_Pion1", "Particle.at(MCtau3pinu_indices[1])")
            .Define("MC_Pion2", "Particle.at(MCtau3pinu_indices[2])")
            .Define("MC_Pion3", "Particle.at(MCtau3pinu_indices[3])")
            ##.Define("MC_tau_CC", "Particle.at(MCtau3pinu_indices_CC[0])")
            .Define("MC_Pion1_CC", "Particle.at(MCtau3pinu_indices_CC[1])")
            .Define("MC_Pion2_CC", "Particle.at(MCtau3pinu_indices_CC[2])")
            .Define("MC_Pion3_CC", "Particle.at(MCtau3pinu_indices_CC[3])")
            .Define("PiMC_M", "VertexAnalysis::tau3pi_MC_mass(MCtau3pinu_indices, Particle)")
            .Define("PiMC_M_CC", "VertexAnalysis::tau3pi_MC_mass(MCtau3pinu_indices_CC, Particle)")
            
            # Decay vertex (an `edm4hep::Vector3d`) of the tau (MC) = production
            .Define("TauMCVertex", "MC_Pion1.vertex")
            .Define("TauMCVertex_x", "TauMCVertex.x")
            .Define("TauMCVertex_y", "TauMCVertex.y")
            .Define("TauMCVertex_z", "TauMCVertex.z")
            .Define("TauMCVertex_CC", "MC_Pion1_CC.vertex")
            .Define("TauMCVertex_x_CC", "TauMCVertex_CC.x")
            .Define("TauMCVertex_y_CC", "TauMCVertex_CC.y")
            .Define("TauMCVertex_z_CC", "TauMCVertex_CC.z")
            
            .Define("ZMCVertex", "MC_tau.vertex")
            .Define("ZMCVertex_x", "ZMCVertex.x")
            .Define("ZMCVertex_y", "ZMCVertex.y")
            .Define("ZMCVertex_z", "ZMCVertex.z")
            ##.Define("dx", "TauMCVertex_x - ZMCVertex_x")
            ##.Define("dy", "TauMCVertex_y - ZMCVertex_y")
            ##.Define("dz", "TauMCVertex_z - ZMCVertex_z")
            ##.Define("dx_CC", "TauMCVertex_x_CC - ZMCVertex_x")
            ##.Define("dy_CC", "TauMCVertex_y_CC - ZMCVertex_y")
            ##.Define("dz_CC", "TauMCVertex_z_CC - ZMCVertex_z")

            .Alias("MCIndex", "MCRecoAssociations#1.index") # Points to MCParticles
            .Alias("RecoIndex", "MCRecoAssociations#0.index") # Points to RecoParticles
            
            ### PHOTONS ###
            
            .Alias("Photon0", "Photon#0.index") # Reconstructed photon indices
            # Filter all MC photons (PDG ID = 22)
            .Define("MC_all_photons", "MCParticle::sel_pdgID(22, false)(Particle)")
            # Filter stable MC photons (genStatus = 1)
            .Define("MC_stable_photons", "MCParticle::sel_genStatus(1)(MC_all_photons)")
            # Vector of energies of each photon in the event
            .Define("MC_photon_e", "MCParticle::get_e(MC_stable_photons)")
            .Define("MC_photon_e_tot", "ROOT::VecOps::Sum(MC_photon_e)")
            .Define("n_MC_photons", "MC_stable_photons.size()")
        )

        ##################################################################################################
        # MC-Kinematics and parents of both: Tau^+ and Tau^-
        ##################################################################################################
        dframe0 = (
            dframe0
            .Define("MC_tau1", "MCParticle::sel_pdgID(15, false)(Particle)")
            .Define("MC_tau2", "MCParticle::sel_pdgID(-15, false)(Particle)")

            .Define("MC_tau_parent", "MCParticle::get_leptons_origin(MC_tau1, Particle, Parents)")
            .Redefine("MC_tau1", "MC_tau1[MC_tau_parent==23]")

            .Define("MC_tau_parent2", "MCParticle::get_leptons_origin(MC_tau2, Particle, Parents)")
            .Redefine("MC_tau2", "MC_tau2[MC_tau_parent2==23]")

            #.Define("MCTauMinus_px","MCParticle::get_px(MC_tau1)[0]")
            #.Define("MCTauMinus_py","MCParticle::get_py(MC_tau1)[0]")
            #.Define("MCTauMinus_pz","MCParticle::get_pz(MC_tau1)[0]")
            .Define("TauMC_E", "MCParticle::get_e(MC_tau1)[0]")

            #.Define("MCTauPlus_px","MCParticle::get_px(MC_tau2)[0]")
            #.Define("MCTauPlus_py","MCParticle::get_py(MC_tau2)[0]")
            #.Define("MCTauPlus_pz","MCParticle::get_pz(MC_tau2)[0]")
            .Define("TauMC_E_CC", "MCParticle::get_e(MC_tau2)[0]")
            
            .Define("PiMC_px", "MC_Pion1.momentum.x + MC_Pion2.momentum.x + MC_Pion3.momentum.x")
            .Define("PiMC_py", "MC_Pion1.momentum.y + MC_Pion2.momentum.y + MC_Pion3.momentum.y")
            .Define("PiMC_pz", "MC_Pion1.momentum.z + MC_Pion2.momentum.z + MC_Pion3.momentum.z")
            .Define("PiMC_px_CC", "MC_Pion1_CC.momentum.x + MC_Pion2_CC.momentum.x + MC_Pion3_CC.momentum.x")
            .Define("PiMC_py_CC", "MC_Pion1_CC.momentum.y + MC_Pion2_CC.momentum.y + MC_Pion3_CC.momentum.y")
            .Define("PiMC_pz_CC", "MC_Pion1_CC.momentum.z + MC_Pion2_CC.momentum.z + MC_Pion3_CC.momentum.z")
        )
        print(f"No. of tau^- -> pi^- pi^+ pi^- nu_tau and CC events (MC): {dframe0.Count().GetValue()}")
        #RDFanalysis.output_var.extend(["MCTauMinus_px","MCTauMinus_py","MCTauMinus_pz","MCTauMinus_E"])
        #RDFanalysis.output_var.extend(["MCTauPlus_px","MCTauPlus_py","MCTauPlus_pz","MCTauPlus_E"])
        RDFanalysis.output_var.extend(["TauMC_E", "TauMC_E_CC"])
        RDFanalysis.output_var.extend(["TauMCVertex_x", "TauMCVertex_y", "TauMCVertex_z"])
        RDFanalysis.output_var.extend(["TauMCVertex_x_CC", "TauMCVertex_y_CC", "TauMCVertex_z_CC"])
        RDFanalysis.output_var.extend(["PiMC_px", "PiMC_py", "PiMC_pz"])
        RDFanalysis.output_var.extend(["PiMC_px_CC", "PiMC_py_CC", "PiMC_pz_CC"])
        RDFanalysis.output_var.extend(["ZMCVertex_x","ZMCVertex_y","ZMCVertex_z"])
        RDFanalysis.output_var.extend(["MC_photon_e_tot", "n_MC_photons"])
        
        ##################################################################################################
        # Reconstructed Kinematics: MC seeded
        ##################################################################################################
        dframe2 = (
            # Returns the RecoParticles associated with the tau decay products.
            # The size of this collection is always 4 provided that
            # MCtau3pinu_indices is not empty, possibly including "dummy"
            # particles in case one of the legs did not make a RecoParticle
            # (e.g. because it is outside the tracker acceptance). This is done
            # on purpose, in order to maintain the mapping with the indices ---
            # i.e. the 1st particle in the list MCtau3pinu_indices is the pi-,
            # , etc.
            # (selRP_matched_to_list ignores the unstable MC particles that are
            # in the input list of indices hence the mother particle, which is
            # the [0] element of the MCtau3pinu_indices vector).
            #
            # The matching between RecoParticles and MCParticles requires 4
            # collections. For more detail, see
            # https://github.com/HEP-FCC/FCCAnalyses/tree/master/examples/basics
            dframe0.Define("TauRecoParticles", "ReconstructedParticle2MC::selRP_matched_to_list(MCtau3pinu_indices, RecoIndex, MCIndex, ReconstructedParticles, Particle)")
            .Define("TauRecoParticles_CC", "ReconstructedParticle2MC::selRP_matched_to_list(MCtau3pinu_indices_CC, RecoIndex, MCIndex, ReconstructedParticles, Particle)")
            
            # the corresponding tracks --- here, dummy particles, if any, are
            # removed, i.e. one may have < 4 tracks, e.g. if one muon or kaon
            # was emitted outside of the acceptance
            .Define("TauTracks", "ReconstructedParticle2Track::getRP2TRK(TauRecoParticles, EFlowTrack_1)")
            .Define("TauTracks_CC", "ReconstructedParticle2Track::getRP2TRK(TauRecoParticles_CC, EFlowTrack_1)")
                    
            # number of tracks in this BsTracks collection (= the #tracks used to reconstruct the Bs vertex)
            .Define("n_TauTracks", "ReconstructedParticle2Track::getTK_n(TauTracks)")
            .Define("n_TauTracks_CC", "ReconstructedParticle2Track::getTK_n(TauTracks_CC)")
                    
            # number of tracks in this BsTracks collection (= the #tracks used to reconstruct the Bs vertex)
            .Filter("n_TauTracks == 3 && n_TauTracks_CC == 3")
            
            # Computes pre-vertex-fit momenta of the 3 pions
            .Define("Tracks_px", "Sum(ReconstructedParticle::get_px(TauRecoParticles))")
            .Define("Tracks_py", "Sum(ReconstructedParticle::get_py(TauRecoParticles))")
            .Define("Tracks_pz", "Sum(ReconstructedParticle::get_pz(TauRecoParticles))")
            .Define("Tracks_px_CC", "Sum(ReconstructedParticle::get_px(TauRecoParticles_CC))")
            .Define("Tracks_py_CC", "Sum(ReconstructedParticle::get_py(TauRecoParticles_CC))")
            .Define("Tracks_pz_CC", "Sum(ReconstructedParticle::get_pz(TauRecoParticles_CC))")
            
            # Fit the tracks to a common vertex. That would be a secondary vertex, hence we put a "2" as the
            # first argument of VertexFitter_Tk: First the full object, of type Vertexing::FCCAnalysesVertex
            .Define("TauVertexObject", "VertexFitterSimple::VertexFitter_Tk(2, TauTracks)")
            .Define("TauVertexObject_CC", "VertexFitterSimple::VertexFitter_Tk(2, TauTracks_CC)")
            # from which we extract the edm4hep::VertexData object, which contains the vertex position in mm

            .Define("TauVertex", "VertexingUtils::get_VertexData(TauVertexObject)")
            .Define("TauVertex_chi2", "TauVertex.chi2")
            .Filter("TauVertex_chi2 < 4.5") # Clean-up for badly reco'ed events
            .Define("TauVertex_x", "TauVertex.position.x")
            .Define("TauVertex_y", "TauVertex.position.y")
            .Define("TauVertex_z", "TauVertex.position.z")
            .Define("TauVertex_CC", "VertexingUtils::get_VertexData(TauVertexObject_CC)")
            .Define("TauVertex_chi2_CC", "TauVertex_CC.chi2")
            .Filter("TauVertex_chi2_CC < 4.5") # Clean-up for badly reco'ed events
            .Define("TauVertex_x_CC", "TauVertex_CC.position.x")
            .Define("TauVertex_y_CC", "TauVertex_CC.position.y")
            .Define("TauVertex_z_CC", "TauVertex_CC.position.z")
            
            .Define("P", "VertexAnalysis::TauVertexMomentum(TauVertexObject)")
            .Define("Pi_px", "P.X()")
            .Define("Pi_py", "P.Y()")
            .Define("Pi_pz", "P.Z()")
            .Define("P_CC", "VertexAnalysis::TauVertexMomentum(TauVertexObject_CC)")
            .Define("Pi_px_CC", "P_CC.X()")
            .Define("Pi_py_CC", "P_CC.Y()")
            .Define("Pi_pz_CC", "P_CC.Z()")
            #.Define("TauVertexCharge","Sum(ReconstructedParticle::get_charge(TauRecoParticles))")
            
            .Define("invMass_3pi_fit", "VertexAnalysis::tau3pi_vertex_mass(TauVertexObject)")
            .Filter("invMass_3pi_fit < 1.8") # Clean-up for badly reco'ed events
            .Define("invMass_3pi_fit_CC", "VertexAnalysis::tau3pi_vertex_mass(TauVertexObject_CC)")
            .Filter("invMass_3pi_fit_CC < 1.8") # Clean-up for badly reco'ed events
            .Define("invMass_3pi_MC", "VertexAnalysis::tau3pi_MC_mass(MCtau3pinu_indices, Particle)")
            .Define("invMass_3pi_MC_CC", "VertexAnalysis::tau3pi_MC_mass(MCtau3pinu_indices_CC, Particle)")
            
            ##.Define("Mass", "VertexAnalysis::TauMassMC(TauMC_E, PionsMC_M, dx, dy, dz, PionsMC_px, PionsMC_py, PionsMC_pz)") 
            ##.Define("Mass_CC", "VertexAnalysis::TauMassMC(TauMC_E_CC, PionsMC_M_CC, dx_CC, dy_CC, dz_CC, PionsMC_px_CC, PionsMC_py_CC, PionsMC_pz_CC)")
            ##.Define("MassMC", "Mass.first")
            ##.Define("MassMC_CC", "Mass_CC.first")
            ##.Define("AngleMC", "mass_cb.second")
            
            # Estimates the primary vertex (PV) based on the reconstructed decay vertices of both taus.
            # 3rd, 4th, 5th arguments: bsc_x, bsc_y, bsc_z in mm, where "bsc" is beam spot constraint (sigma).
            .Define("ZVertex", "VertexAnalysis::EstimatePV(TauVertex_x, TauVertex_y, TauVertex_z, TauVertex_x_CC, TauVertex_y_CC, TauVertex_z_CC, 6e-3, 24e-6, 397e-3)")
            .Define("ZVertex_x", "ZVertex.x")
            .Define("ZVertex_y", "ZVertex.y")
            .Define("ZVertex_z", "ZVertex.z")
            
            # Computes tau mass and cosine of the angle between the tau and the 3-pion system
            .Define("mass_cb_reco", "VertexAnalysis::TauMass(TauVertexObject, TauVertex_x, TauVertex_y, TauVertex_z, ZVertex_x, ZVertex_y, ZVertex_z)")
            .Define("mass_cb_reco_CC", "VertexAnalysis::TauMass(TauVertexObject_CC, TauVertex_x_CC, TauVertex_y_CC, TauVertex_z_CC, ZVertex_x, ZVertex_y, ZVertex_z)")      
            .Define("MassReco", "mass_cb_reco.first")
            .Define("MassReco_CC", "mass_cb_reco_CC.first")
            .Define("AngleReco", "mass_cb_reco.second")
            
            # Computes tau mass using truth values for the 3-pion system momenta
            .Define("mass_cb_reco_truem", "VertexAnalysis::TauMass(PiMC_px, PiMC_py, PiMC_pz, PiMC_M, TauVertex_x, TauVertex_y, TauVertex_z, ZVertex_x, ZVertex_y, ZVertex_z)")
            .Define("MassReco_truem", "mass_cb_reco_truem.first")
            .Define("mass_cb_reco_truem_CC", "VertexAnalysis::TauMass(PiMC_px_CC, PiMC_py_CC, PiMC_pz_CC, PiMC_M_CC, TauVertex_x_CC, TauVertex_y_CC, TauVertex_z_CC, ZVertex_x, ZVertex_y, ZVertex_z)")
            .Define("MassReco_truem_CC", "mass_cb_reco_truem_CC.first")
            
            .Define("TauLifet", "VertexAnalysis::TauLifetime(TauVertex_x, TauVertex_y, TauVertex_z, ZVertex_x, ZVertex_y, ZVertex_z)")
            .Define("TauLifet_CC", "VertexAnalysis::TauLifetime(TauVertex_x_CC, TauVertex_y_CC, TauVertex_z_CC, ZVertex_x, ZVertex_y, ZVertex_z)")
            
            # The "raw" mass --- using the track momenta at their dca:
            ##.Define("invMass_3pi_raw", "VertexAnalysis::tau3pi_raw_mass(TauRecoParticles)")
            ##.Define("invMass_3pi_raw_CC", "VertexAnalysis::tau3pi_raw_mass(TauRecoParticles_CC)")
            ##.Define('RP_MC_index', "ReconstructedParticle2MC::getRP2MC_index(RecoIndex,MCIndex,ReconstructedParticles)")
            
            ### PHOTONS ###
            
            .Define("Photons", "ReconstructedParticle::get(Photon0, ReconstructedParticles)")
            .Define("Photons_E", "ReconstructedParticle::get_e(Photons)")
            .Define("Photons_E_tot", "ROOT::VecOps::Sum(Photons_E)")
            .Define("n_photons", "Photons.size()")
        )
        print(f"No. of reconstructed tau->3pi events (MC-seeded): {dframe2.Count().GetValue()}")
        
        RDFanalysis.output_var.extend(["TauVertex_x", "TauVertex_y", "TauVertex_z"])
        RDFanalysis.output_var.extend(["TauVertex_x_CC", "TauVertex_y_CC", "TauVertex_z_CC"])
        RDFanalysis.output_var.extend(["ZVertex_x", "ZVertex_y", "ZVertex_z"])
        RDFanalysis.output_var.extend(["MassReco", "MassReco_CC"])
        RDFanalysis.output_var.extend(["MassReco_truem", "MassReco_truem_CC"])
        RDFanalysis.output_var.extend(["TauLifet", "TauLifet_CC"])
        RDFanalysis.output_var.extend(["invMass_3pi_MC", "invMass_3pi_MC_CC"])
        #RDFanalysis.output_var.extend(["invMass_3pi_raw", "invMass_3pi_raw_CC"])
        RDFanalysis.output_var.extend(["invMass_3pi_fit", "invMass_3pi_fit_CC"])
        RDFanalysis.output_var.extend(["Pi_px", "Pi_py", "Pi_pz"])
        RDFanalysis.output_var.extend(["Pi_px_CC", "Pi_py_CC", "Pi_pz_CC"])
        RDFanalysis.output_var.extend(["Tracks_px", "Tracks_py", "Tracks_pz"])
        RDFanalysis.output_var.extend(["Tracks_px_CC", "Tracks_py_CC", "Tracks_pz_CC"])
        RDFanalysis.output_var.extend(["Photons_E_tot", "n_photons"])
        
        print("Analysis done :)")
        return dframe2

    def output():
        '''
        Output variables which will be saved to output root file.
        '''

        branch_list=RDFanalysis.output_var

        return branch_list