# REORG_LOG.md

Complete audit log of the reorganization pass, **derived from current filesystem state** (re-hashed at destination) + `_reorg/lib_hashes.json`. No file contents were edited; every action is a move or a byte-faithful copy. **Zero deletions** this pass â€” all removal candidates were *parked* into `archive/_to_delete/` for the user to empty by hand.

Generated: 2026-05-29. Total logged rows: 310.

Columns: action | source path | destination path | MD5-before | MD5-after. For `move` / `park` the source no longer exists at its old path (the operation was a same-volume rename verified byte-faithful inline at action time in stages B & D); before==after is the hash captured now at the destination. For `copy-from-lib` / `copy-as-is` both hashes are computed live (source still present), so the row is a genuine before-vs-after check. `move-bulk-verify` rows carry file-count + total-bytes for the whole directory (binary/data content verified at folder granularity).

## lib/

| action | source | destination | MD5-before | MD5-after |
|---|---|---|---|---|
| copy-to-lib | codes\run_codes\munich\baby_bayes.R | lib\baby_bayes.R | 8EDC5C1BB9E78BAAEA741A7EC2D38C4F | 8EDC5C1BB9E78BAAEA741A7EC2D38C4F |
| copy-to-lib | codes\test and comparison\canonical_schema.R | lib\canonical_schema.R | 76EAB373FE11104E7B7AA7A507D5A104 | 76EAB373FE11104E7B7AA7A507D5A104 |
| copy-to-lib | codes\test and comparison\fit_ours_interaction_wrapper.R | lib\fit_ours_interaction_wrapper.R | DADA086BC3CA29F082533C63F1807A25 | DADA086BC3CA29F082533C63F1807A25 |
| copy-to-lib | codes\run_codes\munich\fit_spatial_reml.R | lib\fit_spatial_reml.R | 0F803D62A7CE9448CF6591977507AD78 | 0F803D62A7CE9448CF6591977507AD78 |
| copy-to-lib | codes\run_codes\munich\gibbs_bayes.R | lib\gibbs_bayes.R | B8D2CD8B3532518C6539F6ABD8EFE841 | B8D2CD8B3532518C6539F6ABD8EFE841 |
| copy-to-lib | codes\run_codes\munich\gibbs_interaction.R | lib\gibbs_interaction.R | 9C0A1288E0E8899AE1ED15FC4BF63AFC | 9C0A1288E0E8899AE1ED15FC4BF63AFC |
| copy-to-lib | codes\run_codes\munich\gibbs_stage_c_full.R | lib\gibbs_stage_c_full.R | D22348BA7CD55F6DC3459BC08E7F387C | D22348BA7CD55F6DC3459BC08E7F387C |
| copy-to-lib | codes\run_codes\munich\ls_basis.R | lib\ls_basis.R | 6A3132B2720451E7C7B53029CE5C0EFE | 6A3132B2720451E7C7B53029CE5C0EFE |
| copy-to-lib | codes\run_codes\munich\ls_interaction.R | lib\ls_interaction.R | 4AE41D9C80024EF4A3732409BE917D47 | 4AE41D9C80024EF4A3732409BE917D47 |
| copy-to-lib | codes\run_codes\munich\ls_interaction_core.cpp | lib\ls_interaction_core.cpp | 81F931C45DE492857D5DE0CE4FEB00E3 | 81F931C45DE492857D5DE0CE4FEB00E3 |
| copy-to-lib | codes\run_codes\munich\marginal_utils.R | lib\marginal_utils.R | 04498862A8B2131139691A94630EB827 | 04498862A8B2131139691A94630EB827 |
| copy-to-lib | codes\run_codes\munich\spatial_utils.R | lib\spatial_utils.R | F0216DB4FD9FD76269EDC0F3E897CD40 | F0216DB4FD9FD76269EDC0F3E897CD40 |

## 01_prototypes

| action | source | destination | MD5-before | MD5-after |
|---|---|---|---|---|
| move | codes\BABY_BAYES_codes\baby_bayes.R | 01_prototypes\baby_bayes\baby_bayes.R | 8EDC5C1BB9E78BAAEA741A7EC2D38C4F | 8EDC5C1BB9E78BAAEA741A7EC2D38C4F |
| move | codes\BABY_BAYES_codes\gibbs_bayes.R | 01_prototypes\baby_bayes\gibbs_bayes.R | B8D2CD8B3532518C6539F6ABD8EFE841 | B8D2CD8B3532518C6539F6ABD8EFE841 |
| move | codes\BABY_BAYES_codes\gibbs_bayes_v2.R | 01_prototypes\baby_bayes\gibbs_bayes_v2.R | 0E90135604812E7B211C5BAEDE64EB26 | 0E90135604812E7B211C5BAEDE64EB26 |
| move | codes\BABY_BAYES_codes\gibbs_stage_c.R | 01_prototypes\baby_bayes\gibbs_stage_c.R | 5727DDA8F7AD49A2DA025439698F32AE | 5727DDA8F7AD49A2DA025439698F32AE |
| move | codes\BABY_BAYES_codes\run_baby_bayes.R | 01_prototypes\baby_bayes\run_baby_bayes.R | 5A87BCEEA67A9F0E93982110985ADB37 | 5A87BCEEA67A9F0E93982110985ADB37 |
| move | codes\BABY_BAYES_codes\run_gibbs_bayes.R | 01_prototypes\baby_bayes\run_gibbs_bayes.R | 54EB19AD372BDF678AA112856DD235E0 | 54EB19AD372BDF678AA112856DD235E0 |
| move | codes\BABY_BAYES_codes\run_gibbs_v2.R | 01_prototypes\baby_bayes\run_gibbs_v2.R | F73C59AEC496BC10172008268007BFAC | F73C59AEC496BC10172008268007BFAC |
| move | codes\BABY_BAYES_codes\run_stage_c.R | 01_prototypes\baby_bayes\run_stage_c.R | 3168528B23A160BCEC6E162CECC57ED7 | 3168528B23A160BCEC6E162CECC57ED7 |
| move | codes\BABY_BAYES_codes\run_stage_c_all.R | 01_prototypes\baby_bayes\run_stage_c_all.R | B30BBE46F1AA6B6C116D6FCFEDAC23B1 | B30BBE46F1AA6B6C116D6FCFEDAC23B1 |
| move | codes\BABY_BAYES_codes\submit_stage_c.sub | 01_prototypes\baby_bayes\submit_stage_c.sub | 0C459E60262CAF317473AF62B6751F5D | 0C459E60262CAF317473AF62B6751F5D |
| move | codes\hist_codes\stage1_ls_only.R | 01_prototypes\hist_codes\stage1_ls_only.R | 5D5CDE33881FFC89A5EB0019BEB81A89 | 5D5CDE33881FFC89A5EB0019BEB81A89 |
| move | codes\hist_codes\stage1b_ls_only.R | 01_prototypes\hist_codes\stage1b_ls_only.R | A3185BFFF6E96540E69CDB67264C56E3 | A3185BFFF6E96540E69CDB67264C56E3 |
| move | codes\hist_codes\stage2_spatialX.R | 01_prototypes\hist_codes\stage2_spatialX.R | 6B0252F9E0D5C56B6154E4DA7327ACE7 | 6B0252F9E0D5C56B6154E4DA7327ACE7 |
| move-bulk-verify | codes\BABY_BAYES_codes | 01_prototypes\baby_bayes | count=10;bytes=94286 | count=10;bytes=94286 |
| move-bulk-verify | codes\hist_codes | 01_prototypes\hist_codes | count=3;bytes=17196 | count=3;bytes=17196 |

## 02_reml

| action | source | destination | MD5-before | MD5-after |
|---|---|---|---|---|
| copy-from-lib | lib\ls_basis.R | 02_reml\ls_basis.R | 6A3132B2720451E7C7B53029CE5C0EFE | 6A3132B2720451E7C7B53029CE5C0EFE |
| copy-from-lib | lib\spatial_utils.R | 02_reml\spatial_utils.R | F0216DB4FD9FD76269EDC0F3E897CD40 | F0216DB4FD9FD76269EDC0F3E897CD40 |
| move | codes\REML_codes\boundary_experiments.R | 02_reml\boundary_experiments.R | 4E6F2CAC217C0C25566DE406F33138C0 | 4E6F2CAC217C0C25566DE406F33138C0 |
| move | codes\REML_codes\fit_spatial_reml.R | 02_reml\fit_spatial_reml.R | 0F803D62A7CE9448CF6591977507AD78 | 0F803D62A7CE9448CF6591977507AD78 |
| move | codes\REML_codes\marginal_utils.R | 02_reml\marginal_utils.R | 04498862A8B2131139691A94630EB827 | 04498862A8B2131139691A94630EB827 |
| move | codes\REML_codes\plot_utils.R | 02_reml\plot_utils.R | 608D037E5AC8FE2651DA1C5B3FD69398 | 608D037E5AC8FE2651DA1C5B3FD69398 |
| move | codes\REML_codes\run.R | 02_reml\run.R | 5B569E59C86FF2627DBD81E068418906 | 5B569E59C86FF2627DBD81E068418906 |
| move | codes\REML_codes\run_boundary.R | 02_reml\run_boundary.R | 9A82E1B1A4EF4B77388F5BE0766D9D0F | 9A82E1B1A4EF4B77388F5BE0766D9D0F |
| move | codes\REML_codes\run_boundary_exp1.R | 02_reml\run_boundary_exp1.R | 9B6592602841D3009525190A1CAD6B73 | 9B6592602841D3009525190A1CAD6B73 |
| move | codes\REML_codes\run_sim1_marginal.R | 02_reml\run_sim1_marginal.R | D0B34BABAFAECEAB6E885DC28C6C0146 | D0B34BABAFAECEAB6E885DC28C6C0146 |
| move | codes\REML_codes\run_sim2_marginal.R | 02_reml\run_sim2_marginal.R | 72708EC03A4872E8E12424B57F4BEAAA | 72708EC03A4872E8E12424B57F4BEAAA |
| move | codes\REML_codes\run_sim3_marginal.R | 02_reml\run_sim3_marginal.R | 1CE266E09EE4C8367A1DD0459695E100 | 1CE266E09EE4C8367A1DD0459695E100 |
| move | codes\REML_codes\run_sim3_marginal_n1000.R | 02_reml\run_sim3_marginal_n1000.R | 863AD74DD5F770BE23FB58DDF2FCFCFE | 863AD74DD5F770BE23FB58DDF2FCFCFE |
| move | codes\REML_codes\sim1_sweep.R | 02_reml\sim1_sweep.R | 9129F3B715646575CD22007FBA50F2C6 | 9129F3B715646575CD22007FBA50F2C6 |
| move | codes\REML_codes\sim2_sweep.R | 02_reml\sim2_sweep.R | ADD67BBF891A36C39513607E2D4B4BB9 | ADD67BBF891A36C39513607E2D4B4BB9 |
| move | codes\REML_codes\Sim3_plot_surfaces.R | 02_reml\Sim3_plot_surfaces.R | 8C11EA81B1B0BC4C575590E8A4D9DF0B | 8C11EA81B1B0BC4C575590E8A4D9DF0B |
| move | codes\REML_codes\sim3_sweep.R | 02_reml\sim3_sweep.R | 4BA02B733A027ED24696063222328F5B | 4BA02B733A027ED24696063222328F5B |
| move | codes\REML_codes\sim4_sweep.R | 02_reml\sim4_sweep.R | 023AD60E642E7102B77AA24AD45691CA | 023AD60E642E7102B77AA24AD45691CA |
| move | codes\REML_codes\sim5_sweep.R | 02_reml\sim5_sweep.R | 494AB74A0DDA06FCA618FD8EA56BC273 | 494AB74A0DDA06FCA618FD8EA56BC273 |
| move | submit_boundary.sub | 02_reml\submit_boundary.sub | 93B43DBCC62CCC54A4440972198DF987 | 93B43DBCC62CCC54A4440972198DF987 |
| move | submit_exp1.sub | 02_reml\submit_exp1.sub | 3086305FAEBE9C85F539426E04C8AE66 | 3086305FAEBE9C85F539426E04C8AE66 |
| move | codes\REML_codes\variable_selection.R | 02_reml\variable_selection.R | 1B452A3418C9FC4ED44B20F5935A8168 | 1B452A3418C9FC4ED44B20F5935A8168 |
| move-bulk-verify | codes\REML_codes | 02_reml | count=22;bytes=121911 | count=22;bytes=121911 |

## 03_bayes_main

| action | source | destination | MD5-before | MD5-after |
|---|---|---|---|---|
| copy-from-lib | lib\baby_bayes.R | 03_bayes_main\n1000\baby_bayes.R | 8EDC5C1BB9E78BAAEA741A7EC2D38C4F | 8EDC5C1BB9E78BAAEA741A7EC2D38C4F |
| copy-from-lib | lib\fit_spatial_reml.R | 03_bayes_main\n1000\fit_spatial_reml.R | 0F803D62A7CE9448CF6591977507AD78 | 0F803D62A7CE9448CF6591977507AD78 |
| copy-from-lib | lib\gibbs_bayes.R | 03_bayes_main\n1000\gibbs_bayes.R | B8D2CD8B3532518C6539F6ABD8EFE841 | B8D2CD8B3532518C6539F6ABD8EFE841 |
| copy-from-lib | lib\gibbs_stage_c_full.R | 03_bayes_main\n1000\gibbs_stage_c_full.R | D22348BA7CD55F6DC3459BC08E7F387C | D22348BA7CD55F6DC3459BC08E7F387C |
| copy-from-lib | lib\ls_basis.R | 03_bayes_main\n1000\ls_basis.R | 6A3132B2720451E7C7B53029CE5C0EFE | 6A3132B2720451E7C7B53029CE5C0EFE |
| copy-from-lib | lib\marginal_utils.R | 03_bayes_main\n1000\marginal_utils.R | 04498862A8B2131139691A94630EB827 | 04498862A8B2131139691A94630EB827 |
| copy-from-lib | lib\spatial_utils.R | 03_bayes_main\n1000\spatial_utils.R | F0216DB4FD9FD76269EDC0F3E897CD40 | F0216DB4FD9FD76269EDC0F3E897CD40 |
| copy-from-lib | lib\baby_bayes.R | 03_bayes_main\n400\baby_bayes.R | 8EDC5C1BB9E78BAAEA741A7EC2D38C4F | 8EDC5C1BB9E78BAAEA741A7EC2D38C4F |
| copy-from-lib | lib\fit_spatial_reml.R | 03_bayes_main\n400\fit_spatial_reml.R | 0F803D62A7CE9448CF6591977507AD78 | 0F803D62A7CE9448CF6591977507AD78 |
| copy-from-lib | lib\gibbs_bayes.R | 03_bayes_main\n400\gibbs_bayes.R | B8D2CD8B3532518C6539F6ABD8EFE841 | B8D2CD8B3532518C6539F6ABD8EFE841 |
| copy-from-lib | lib\ls_basis.R | 03_bayes_main\n400\ls_basis.R | 6A3132B2720451E7C7B53029CE5C0EFE | 6A3132B2720451E7C7B53029CE5C0EFE |
| copy-from-lib | lib\marginal_utils.R | 03_bayes_main\n400\marginal_utils.R | 04498862A8B2131139691A94630EB827 | 04498862A8B2131139691A94630EB827 |
| copy-from-lib | lib\spatial_utils.R | 03_bayes_main\n400\spatial_utils.R | F0216DB4FD9FD76269EDC0F3E897CD40 | F0216DB4FD9FD76269EDC0F3E897CD40 |
| move | codes\bayes_n1000_mcmc10000_codes\bayes_overlay\run_sim1_n1000_overlay.R | 03_bayes_main\n1000\bayes_overlay\run_sim1_n1000_overlay.R | A7A18AC6F6AA8CCB195CF0F599904B02 | A7A18AC6F6AA8CCB195CF0F599904B02 |
| move | codes\bayes_n1000_mcmc10000_codes\bayes_overlay\run_sim2_n1000_overlay.R | 03_bayes_main\n1000\bayes_overlay\run_sim2_n1000_overlay.R | 4550E0B353215DC6BFCE27C58F5BA2FD | 4550E0B353215DC6BFCE27C58F5BA2FD |
| move | codes\bayes_n1000_mcmc10000_codes\bayes_overlay\run_sim3_n1000_overlay.R | 03_bayes_main\n1000\bayes_overlay\run_sim3_n1000_overlay.R | F243845C6C4C3763653C446194F0A8D0 | F243845C6C4C3763653C446194F0A8D0 |
| move | codes\bayes_n1000_mcmc10000_codes\bayes_overlay\run_sim4_n1000_overlay.R | 03_bayes_main\n1000\bayes_overlay\run_sim4_n1000_overlay.R | DF1FD35E2A3A9977451AE410061CA41C | DF1FD35E2A3A9977451AE410061CA41C |
| move | codes\bayes_n1000_mcmc10000_codes\bayes_overlay\run_sim4_plus_n1000_overlay.R | 03_bayes_main\n1000\bayes_overlay\run_sim4_plus_n1000_overlay.R | 1CD19BD8D5034A5CF4E918B825181952 | 1CD19BD8D5034A5CF4E918B825181952 |
| move | codes\bayes_n1000_mcmc10000_codes\bayes_overlay\sim2_n1000_overlay.sub | 03_bayes_main\n1000\bayes_overlay\sim2_n1000_overlay.sub | CB797E8726F546639C20D2839D33FEB7 | CB797E8726F546639C20D2839D33FEB7 |
| move | codes\bayes_n1000_mcmc10000_codes\bayes_overlay\sim3_n1000_overlay.sub | 03_bayes_main\n1000\bayes_overlay\sim3_n1000_overlay.sub | 8803B31585E421758080D7D95466DF24 | 8803B31585E421758080D7D95466DF24 |
| move | codes\bayes_n1000_mcmc10000_codes\bayes_overlay\sim4_n1000_overlay.sub | 03_bayes_main\n1000\bayes_overlay\sim4_n1000_overlay.sub | 2374BAF41416DAC6DB83B0C5F7D57B94 | 2374BAF41416DAC6DB83B0C5F7D57B94 |
| move | codes\bayes_n1000_mcmc10000_codes\bayes_overlay\sim4plus_n1000_overlay.sub | 03_bayes_main\n1000\bayes_overlay\sim4plus_n1000_overlay.sub | 6A458A6BA2749B118D9C14B6FB9D9B18 | 6A458A6BA2749B118D9C14B6FB9D9B18 |
| move | codes\bayes_n1000_mcmc10000_codes\bayes_test_sim1\run_sim1_n1000_plot_new.R | 03_bayes_main\n1000\bayes_test_sim1\run_sim1_n1000_plot_new.R | 334FFB0D8252C2F9B512DF1B47CBB735 | 334FFB0D8252C2F9B512DF1B47CBB735 |
| move | codes\bayes_n1000_mcmc10000_codes\bayes_test_sim1\submit_sim1_n1000_overlay.sub | 03_bayes_main\n1000\bayes_test_sim1\submit_sim1_n1000_overlay.sub | 8389693A6523D30FFD402EA69A8C8F7C | 8389693A6523D30FFD402EA69A8C8F7C |
| move | codes\bayes_n1000_mcmc10000_codes\bayes_test_sim1\submit_sim1_n1000_plot_new.sub | 03_bayes_main\n1000\bayes_test_sim1\submit_sim1_n1000_plot_new.sub | C3949B3F5900B602389DFD686E1F446F | C3949B3F5900B602389DFD686E1F446F |
| move | codes\bayes_n1000_mcmc10000_codes\run_bayes_sim4_n1000.R | 03_bayes_main\n1000\run_bayes_sim4_n1000.R | 0B7E3140F083CA968DD9423EAD06C59E | 0B7E3140F083CA968DD9423EAD06C59E |
| move | codes\bayes_n1000_mcmc10000_codes\run_prior_sensitivity_sim4.R | 03_bayes_main\n1000\run_prior_sensitivity_sim4.R | 9C6D8A5984E7ED4B60D4F22DC74B4FBC | 9C6D8A5984E7ED4B60D4F22DC74B4FBC |
| move | codes\bayes_n1000_mcmc10000_codes\run_sample_size_comparison.R | 03_bayes_main\n1000\run_sample_size_comparison.R | E5F2D7EC52BF1236E8A9D954E817869E | E5F2D7EC52BF1236E8A9D954E817869E |
| move | codes\bayes_n1000_mcmc10000_codes\run_sim1_n1000.R | 03_bayes_main\n1000\run_sim1_n1000.R | 292CF7463F5CE164741C4044A695941E | 292CF7463F5CE164741C4044A695941E |
| move | codes\bayes_n1000_mcmc10000_codes\run_sim2_n1000.R | 03_bayes_main\n1000\run_sim2_n1000.R | 8AF0A260934091E72B98E46DAF7412DC | 8AF0A260934091E72B98E46DAF7412DC |
| move | codes\bayes_n1000_mcmc10000_codes\run_sim3_n1000.R | 03_bayes_main\n1000\run_sim3_n1000.R | 74D57EFDAB84A147DE158F5BAF747EB4 | 74D57EFDAB84A147DE158F5BAF747EB4 |
| move | codes\bayes_n1000_mcmc10000_codes\run_sim4_n1000.R | 03_bayes_main\n1000\run_sim4_n1000.R | 95207B402BBFD9BEEC3ADD7BE3C3BA9A | 95207B402BBFD9BEEC3ADD7BE3C3BA9A |
| move | codes\bayes_n1000_mcmc10000_codes\run_sim4_plus_n1000.R | 03_bayes_main\n1000\run_sim4_plus_n1000.R | 4DB647EAE2796EF8C434A02F21A94732 | 4DB647EAE2796EF8C434A02F21A94732 |
| move | codes\bayes_n1000_mcmc10000_codes\submit_bayes_sim4_n1000.sub | 03_bayes_main\n1000\submit_bayes_sim4_n1000.sub | 2F7700E7C4622A507261C736FC63DEF5 | 2F7700E7C4622A507261C736FC63DEF5 |
| move | codes\bayes_n1000_mcmc10000_codes\submit_prior_sensitivity.sub | 03_bayes_main\n1000\submit_prior_sensitivity.sub | 5ACAAA51F09A6A6374E57842C8A30FD0 | 5ACAAA51F09A6A6374E57842C8A30FD0 |
| move | codes\bayes_n1000_mcmc10000_codes\submit_sample_size_comparison.sub | 03_bayes_main\n1000\submit_sample_size_comparison.sub | 87347C09A9367784838B08470329C521 | 87347C09A9367784838B08470329C521 |
| move | codes\bayes_n1000_mcmc10000_codes\submit_sim1_n1000.sub | 03_bayes_main\n1000\submit_sim1_n1000.sub | 77887C2259694C92F1D08B0B808D455D | 77887C2259694C92F1D08B0B808D455D |
| move | codes\bayes_n1000_mcmc10000_codes\submit_sim2_n1000.sub | 03_bayes_main\n1000\submit_sim2_n1000.sub | DBF226A2902CA111E2A360A6A0ACA4EC | DBF226A2902CA111E2A360A6A0ACA4EC |
| move | codes\bayes_n1000_mcmc10000_codes\submit_sim3_n1000.sub | 03_bayes_main\n1000\submit_sim3_n1000.sub | 809A65FF20AAE0210531D3D45A2E18F8 | 809A65FF20AAE0210531D3D45A2E18F8 |
| move | codes\bayes_n1000_mcmc10000_codes\submit_sim4_n1000_MH1.sub | 03_bayes_main\n1000\submit_sim4_n1000_MH1.sub | ABBA809C8DFA4F1EE1EBE1E9AEC30BBB | ABBA809C8DFA4F1EE1EBE1E9AEC30BBB |
| move | codes\bayes_n1000_mcmc10000_codes\submit_sim4_plus_n1000.sub | 03_bayes_main\n1000\submit_sim4_plus_n1000.sub | D562959A2099FCD0D55CE044F36ED8E8 | D562959A2099FCD0D55CE044F36ED8E8 |
| move | codes\bayes_n400_codes\gibbs_stage_c_full.R | 03_bayes_main\n400\gibbs_stage_c_full.R | D22348BA7CD55F6DC3459BC08E7F387C | D22348BA7CD55F6DC3459BC08E7F387C |
| move | codes\bayes_n400_codes\run_full_bayes_all.R | 03_bayes_main\n400\run_full_bayes_all.R | F3B2E1AB0CDDC653FCC428418B2816B5 | F3B2E1AB0CDDC653FCC428418B2816B5 |
| move | codes\bayes_n400_codes\run_multi_seed.R | 03_bayes_main\n400\run_multi_seed.R | F377108CFFC4B32A433FEC5D555280E1 | F377108CFFC4B32A433FEC5D555280E1 |
| move | codes\bayes_n400_codes\run_prior_sensitivity.R | 03_bayes_main\n400\run_prior_sensitivity.R | FA95B86E952CC458DE31CB8E1CCBC712 | FA95B86E952CC458DE31CB8E1CCBC712 |
| move | codes\bayes_n400_codes\run_reml_vs_bayes.R | 03_bayes_main\n400\run_reml_vs_bayes.R | 405FB752753CB7F592863A4C151FEA70 | 405FB752753CB7F592863A4C151FEA70 |
| move | codes\bayes_n400_codes\submit_full_bayes.sub | 03_bayes_main\n400\submit_full_bayes.sub | 108DA00515CAD7197EBD7AEA5859C59B | 108DA00515CAD7197EBD7AEA5859C59B |
| move | codes\bayes_n400_codes\submit_multi_seed.sub | 03_bayes_main\n400\submit_multi_seed.sub | CEC3F3EA0F7BAA45CC79E0FF44041942 | CEC3F3EA0F7BAA45CC79E0FF44041942 |
| move | codes\bayes_n400_codes\submit_prior_sensitivity.sub | 03_bayes_main\n400\submit_prior_sensitivity.sub | 5F9E404B3558499D0A6CDD353E8C127D | 5F9E404B3558499D0A6CDD353E8C127D |
| move | codes\bayes_n400_codes\submit_reml_vs_bayes.sub | 03_bayes_main\n400\submit_reml_vs_bayes.sub | 0899FC0D8C19B6CA40EEE23EC4724458 | 0899FC0D8C19B6CA40EEE23EC4724458 |
| move-bulk-verify | codes\bayes_n1000_mcmc10000_codes | 03_bayes_main\n1000 | count=59;bytes=835594 | count=59;bytes=835594 |
| move-bulk-verify | codes\bayes_n400_codes | 03_bayes_main\n400 | count=15;bytes=97071 | count=15;bytes=97071 |

## 04_interaction

| action | source | destination | MD5-before | MD5-after |
|---|---|---|---|---|
| copy-as-is | archive\interaction_old\comparison_stage1\scenario_A\compute_recovery_metrics.R | 04_interaction\scenarioA_compare\compute_recovery_metrics.R | 859E82C58A8BFC034922F344A7B4F0D6 | 859E82C58A8BFC034922F344A7B4F0D6 |
| copy-as-is | archive\interaction_old\comparison_stage1\scenario_A\fit_inla_wrapper.R | 04_interaction\scenarioA_compare\fit_inla_wrapper.R | 2B0B8A4DBE48A2999BEDA7C8B20D7FD5 | 2B0B8A4DBE48A2999BEDA7C8B20D7FD5 |
| copy-as-is | archive\interaction_old\comparison_stage1\scenario_A\fit_mgcv_wrapper.R | 04_interaction\scenarioA_compare\fit_mgcv_wrapper.R | C05F46CA60878A2F272255F5B13D9787 | C05F46CA60878A2F272255F5B13D9787 |
| copy-as-is | archive\interaction_old\comparison_stage1\scenario_A\fit_ours_wrapper.R | 04_interaction\scenarioA_compare\fit_ours_wrapper.R | DB67A88AA37A481DD696F04905F293A7 | DB67A88AA37A481DD696F04905F293A7 |
| copy-as-is | archive\interaction_old\comparison_stage1\scenario_A\run_stage1_scenarioA_multiseed.R | 04_interaction\scenarioA_compare\run_stage1_scenarioA_multiseed.R | 51E212C879591ECF140EBB1A22A5381C | 51E212C879591ECF140EBB1A22A5381C |
| copy-as-is | archive\interaction_old\comparison_stage1\scenario_A\scenarioA_settings.md | 04_interaction\scenarioA_compare\scenarioA_settings.md | F9F885B826B1FAE8581FBAC3446FB19C | F9F885B826B1FAE8581FBAC3446FB19C |
| copy-as-is | archive\interaction_old\comparison_stage1\scenario_A\simulate_scenario_A.R | 04_interaction\scenarioA_compare\simulate_scenario_A.R | BB45788287D6CE2272AEFA40D3706E80 | BB45788287D6CE2272AEFA40D3706E80 |
| copy-as-is | archive\interaction_old\comparison_stage1\scenario_A\test_fit_inla.R | 04_interaction\scenarioA_compare\test_fit_inla.R | 2A027C9C67531F87ECD86384768D3426 | 2A027C9C67531F87ECD86384768D3426 |
| copy-as-is | archive\interaction_old\comparison_stage1\scenario_A\test_fit_mgcv.R | 04_interaction\scenarioA_compare\test_fit_mgcv.R | 53D3FF54D8F5C0FDF487935FA1C4119D | 53D3FF54D8F5C0FDF487935FA1C4119D |
| copy-as-is | archive\interaction_old\comparison_stage1\scenario_A\test_fit_ours.R | 04_interaction\scenarioA_compare\test_fit_ours.R | 833971C6EBC87DBA1A3754832A6103D7 | 833971C6EBC87DBA1A3754832A6103D7 |
| copy-as-is | archive\interaction_old\comparison_stage1\scenario_A\test_simulate_scenario_A.R | 04_interaction\scenarioA_compare\test_simulate_scenario_A.R | 76731691C0A6380C4F42547C7D6B4434 | 76731691C0A6380C4F42547C7D6B4434 |
| copy-from-lib | lib\canonical_schema.R | 04_interaction\pipeline\canonical_schema.R | 76EAB373FE11104E7B7AA7A507D5A104 | 76EAB373FE11104E7B7AA7A507D5A104 |
| copy-from-lib | lib\fit_ours_interaction_wrapper.R | 04_interaction\pipeline\fit_ours_interaction_wrapper.R | DADA086BC3CA29F082533C63F1807A25 | DADA086BC3CA29F082533C63F1807A25 |
| copy-from-lib | lib\fit_spatial_reml.R | 04_interaction\pipeline\fit_spatial_reml.R | 0F803D62A7CE9448CF6591977507AD78 | 0F803D62A7CE9448CF6591977507AD78 |
| copy-from-lib | lib\gibbs_interaction.R | 04_interaction\pipeline\gibbs_interaction.R | 9C0A1288E0E8899AE1ED15FC4BF63AFC | 9C0A1288E0E8899AE1ED15FC4BF63AFC |
| copy-from-lib | lib\gibbs_stage_c_full.R | 04_interaction\pipeline\gibbs_stage_c_full.R | D22348BA7CD55F6DC3459BC08E7F387C | D22348BA7CD55F6DC3459BC08E7F387C |
| copy-from-lib | lib\ls_basis.R | 04_interaction\pipeline\ls_basis.R | 6A3132B2720451E7C7B53029CE5C0EFE | 6A3132B2720451E7C7B53029CE5C0EFE |
| copy-from-lib | lib\ls_interaction.R | 04_interaction\pipeline\ls_interaction.R | 4AE41D9C80024EF4A3732409BE917D47 | 4AE41D9C80024EF4A3732409BE917D47 |
| copy-from-lib | lib\ls_interaction_core.cpp | 04_interaction\pipeline\ls_interaction_core.cpp | 81F931C45DE492857D5DE0CE4FEB00E3 | 81F931C45DE492857D5DE0CE4FEB00E3 |
| copy-from-lib | lib\marginal_utils.R | 04_interaction\pipeline\marginal_utils.R | 04498862A8B2131139691A94630EB827 | 04498862A8B2131139691A94630EB827 |
| copy-from-lib | lib\spatial_utils.R | 04_interaction\pipeline\spatial_utils.R | F0216DB4FD9FD76269EDC0F3E897CD40 | F0216DB4FD9FD76269EDC0F3E897CD40 |
| copy-from-lib | lib\canonical_schema.R | 04_interaction\scenarioA_compare\canonical_schema.R | 76EAB373FE11104E7B7AA7A507D5A104 | 76EAB373FE11104E7B7AA7A507D5A104 |
| copy-from-lib | lib\fit_spatial_reml.R | 04_interaction\scenarioA_compare\fit_spatial_reml.R | 0F803D62A7CE9448CF6591977507AD78 | 0F803D62A7CE9448CF6591977507AD78 |
| copy-from-lib | lib\gibbs_stage_c_full.R | 04_interaction\scenarioA_compare\gibbs_stage_c_full.R | D22348BA7CD55F6DC3459BC08E7F387C | D22348BA7CD55F6DC3459BC08E7F387C |
| copy-from-lib | lib\ls_basis.R | 04_interaction\scenarioA_compare\ls_basis.R | 6A3132B2720451E7C7B53029CE5C0EFE | 6A3132B2720451E7C7B53029CE5C0EFE |
| copy-from-lib | lib\marginal_utils.R | 04_interaction\scenarioA_compare\marginal_utils.R | 04498862A8B2131139691A94630EB827 | 04498862A8B2131139691A94630EB827 |
| copy-from-lib | lib\spatial_utils.R | 04_interaction\scenarioA_compare\spatial_utils.R | F0216DB4FD9FD76269EDC0F3E897CD40 | F0216DB4FD9FD76269EDC0F3E897CD40 |
| move | codes\test and comparison\canonical_schema.R | 04_interaction\scenarioB_tests\canonical_schema.R | 76EAB373FE11104E7B7AA7A507D5A104 | 76EAB373FE11104E7B7AA7A507D5A104 |
| move | codes\test and comparison\compare_h3_metrics.R | 04_interaction\scenarioB_tests\compare_h3_metrics.R | 265AB2E772DCC100E351A9D65812631C | 265AB2E772DCC100E351A9D65812631C |
| move | codes\test and comparison\compare_postfix_multiseed.R | 04_interaction\scenarioB_tests\compare_postfix_multiseed.R | E8F1628F779C99578ED31CA56FEEC3F8 | E8F1628F779C99578ED31CA56FEEC3F8 |
| move | codes\test and comparison\compare_seed1_metrics.R | 04_interaction\scenarioB_tests\compare_seed1_metrics.R | 54D6F3F0FAD37737881C51DB3C031B46 | 54D6F3F0FAD37737881C51DB3C031B46 |
| move | codes\test and comparison\diagnose_orthogonalization.R | 04_interaction\scenarioB_tests\diagnose_orthogonalization.R | D7590EB6C66A0AC285E110DACDC3687D | D7590EB6C66A0AC285E110DACDC3687D |
| move | codes\test and comparison\fit_ours_interaction_wrapper.R | 04_interaction\scenarioB_tests\fit_ours_interaction_wrapper.R | DADA086BC3CA29F082533C63F1807A25 | DADA086BC3CA29F082533C63F1807A25 |
| move | codes\test and comparison\fit_spatial_reml.R | 04_interaction\scenarioB_tests\fit_spatial_reml.R | 0F803D62A7CE9448CF6591977507AD78 | 0F803D62A7CE9448CF6591977507AD78 |
| move | codes\test and comparison\gibbs_interaction.R | 04_interaction\scenarioB_tests\gibbs_interaction.R | 9C0A1288E0E8899AE1ED15FC4BF63AFC | 9C0A1288E0E8899AE1ED15FC4BF63AFC |
| move | codes\test and comparison\gibbs_stage_c_full.R | 04_interaction\scenarioB_tests\gibbs_stage_c_full.R | D22348BA7CD55F6DC3459BC08E7F387C | D22348BA7CD55F6DC3459BC08E7F387C |
| move | codes\test and comparison\ls_basis.R | 04_interaction\scenarioB_tests\ls_basis.R | 6A3132B2720451E7C7B53029CE5C0EFE | 6A3132B2720451E7C7B53029CE5C0EFE |
| move | codes\test and comparison\ls_interaction.R | 04_interaction\scenarioB_tests\ls_interaction.R | 4AE41D9C80024EF4A3732409BE917D47 | 4AE41D9C80024EF4A3732409BE917D47 |
| move | codes\test and comparison\ls_interaction_core.cpp | 04_interaction\scenarioB_tests\ls_interaction_core.cpp | 81F931C45DE492857D5DE0CE4FEB00E3 | 81F931C45DE492857D5DE0CE4FEB00E3 |
| move | codes\test and comparison\meeting_comparison.R | 04_interaction\scenarioB_tests\meeting_comparison.R | 96A71053989D7262FFA1E183E99EDBC0 | 96A71053989D7262FFA1E183E99EDBC0 |
| move | codes\test and comparison\plot_seed1_comparison.R | 04_interaction\scenarioB_tests\plot_seed1_comparison.R | F6DF6E3419B5685E4C20507680C63144 | F6DF6E3419B5685E4C20507680C63144 |
| move | codes\test and comparison\run_mgcv_scenarioB_multiseed.R | 04_interaction\scenarioB_tests\run_mgcv_scenarioB_multiseed.R | 9F3FFD1385F81A70B0C29849AFB58D54 | 9F3FFD1385F81A70B0C29849AFB58D54 |
| move | codes\test and comparison\run_scenarioB_clamp_all_h3.R | 04_interaction\scenarioB_tests\run_scenarioB_clamp_all_h3.R | CCF28B0A8800B1BDF7A6A19888AAE08D | CCF28B0A8800B1BDF7A6A19888AAE08D |
| move | codes\test and comparison\run_scenarioB_clamp_full_h3.R | 04_interaction\scenarioB_tests\run_scenarioB_clamp_full_h3.R | 2922C358129B1CC25CA80548B5C3A99E | 2922C358129B1CC25CA80548B5C3A99E |
| move | codes\test and comparison\run_scenarioB_clamp_main_int_h3.R | 04_interaction\scenarioB_tests\run_scenarioB_clamp_main_int_h3.R | B8E47170AAB95358D02B9F0803079C61 | B8E47170AAB95358D02B9F0803079C61 |
| move | codes\test and comparison\run_scenarioB_confirm_h3.R | 04_interaction\scenarioB_tests\run_scenarioB_confirm_h3.R | 5E9FA7EF55CF8AED41733F741DD54422 | 5E9FA7EF55CF8AED41733F741DD54422 |
| move | codes\test and comparison\run_scenarioB_postfix.R | 04_interaction\scenarioB_tests\run_scenarioB_postfix.R | 712C7D4B40DC4B79D49FF59211B442EA | 712C7D4B40DC4B79D49FF59211B442EA |
| move | codes\test and comparison\run_scenarioB_postfix_multiseed.R | 04_interaction\scenarioB_tests\run_scenarioB_postfix_multiseed.R | 29905560B219D516F2F435A738753E1B | 29905560B219D516F2F435A738753E1B |
| move | codes\test and comparison\run_scenarioB_postfix_orthoFALSE.R | 04_interaction\scenarioB_tests\run_scenarioB_postfix_orthoFALSE.R | E91E38849989063CA69045FEEB6FA433 | E91E38849989063CA69045FEEB6FA433 |
| move | codes\test and comparison\run_scenarioB_tightbs.R | 04_interaction\scenarioB_tests\run_scenarioB_tightbs.R | B452061105BBF48A798285D967EB4074 | B452061105BBF48A798285D967EB4074 |
| move | codes\test and comparison\simulate_scenario_A.R | 04_interaction\scenarioB_tests\simulate_scenario_A.R | BB45788287D6CE2272AEFA40D3706E80 | BB45788287D6CE2272AEFA40D3706E80 |
| move | codes\test and comparison\simulate_scenario_B.R | 04_interaction\scenarioB_tests\simulate_scenario_B.R | 45BE1ACA0C8431684B8E442726D582CC | 45BE1ACA0C8431684B8E442726D582CC |
| move | codes\test and comparison\spatial_utils.R | 04_interaction\scenarioB_tests\spatial_utils.R | F0216DB4FD9FD76269EDC0F3E897CD40 | F0216DB4FD9FD76269EDC0F3E897CD40 |
| move | codes\test and comparison\test_fit_ours_interaction.R | 04_interaction\scenarioB_tests\test_fit_ours_interaction.R | 56A2774E145D3D39A73DD73E5928BFC6 | 56A2774E145D3D39A73DD73E5928BFC6 |
| move | codes\test and comparison\verify_offby1_fix.R | 04_interaction\scenarioB_tests\verify_offby1_fix.R | 36C0A854570FEB945417E200EF03FDFA | 36C0A854570FEB945417E200EF03FDFA |
| move-bulk-verify | codes\test and comparison | 04_interaction\scenarioB_tests | count=74;bytes=1062047858 | count=74;bytes=1062047858 |

## 05_munich_rent

| action | source | destination | MD5-before | MD5-after |
|---|---|---|---|---|
| copy-from-lib | lib\baby_bayes.R | 05_munich_rent\hpc\codes\baby_bayes.R | 8EDC5C1BB9E78BAAEA741A7EC2D38C4F | 8EDC5C1BB9E78BAAEA741A7EC2D38C4F |
| copy-from-lib | lib\fit_spatial_reml.R | 05_munich_rent\hpc\codes\fit_spatial_reml.R | 0F803D62A7CE9448CF6591977507AD78 | 0F803D62A7CE9448CF6591977507AD78 |
| copy-from-lib | lib\gibbs_bayes.R | 05_munich_rent\hpc\codes\gibbs_bayes.R | B8D2CD8B3532518C6539F6ABD8EFE841 | B8D2CD8B3532518C6539F6ABD8EFE841 |
| copy-from-lib | lib\gibbs_stage_c_full.R | 05_munich_rent\hpc\codes\gibbs_stage_c_full.R | D22348BA7CD55F6DC3459BC08E7F387C | D22348BA7CD55F6DC3459BC08E7F387C |
| copy-from-lib | lib\ls_basis.R | 05_munich_rent\hpc\codes\ls_basis.R | 6A3132B2720451E7C7B53029CE5C0EFE | 6A3132B2720451E7C7B53029CE5C0EFE |
| copy-from-lib | lib\marginal_utils.R | 05_munich_rent\hpc\codes\marginal_utils.R | 04498862A8B2131139691A94630EB827 | 04498862A8B2131139691A94630EB827 |
| copy-from-lib | lib\spatial_utils.R | 05_munich_rent\hpc\codes\spatial_utils.R | F0216DB4FD9FD76269EDC0F3E897CD40 | F0216DB4FD9FD76269EDC0F3E897CD40 |
| move | codes\realtest_code\munich_rent\codes\run_munich_rent.R | 05_munich_rent\hpc\codes\run_munich_rent.R | 75FB1165E6A212D9634E320FED12D657 | 75FB1165E6A212D9634E320FED12D657 |
| move | codes\realtest_code\munich_rent\codes\run_munich_rent_n1000.R | 05_munich_rent\hpc\codes\run_munich_rent_n1000.R | F4F3D58EAE6C8357DA6AF171AFBBC308 | F4F3D58EAE6C8357DA6AF171AFBBC308 |
| move | codes\realtest_code\munich_rent\codes\run_munich_rent_n1000_iter5000_testMH.R | 05_munich_rent\hpc\codes\run_munich_rent_n1000_iter5000_testMH.R | B5357EDC43C201897030C4E79B6EF1E1 | B5357EDC43C201897030C4E79B6EF1E1 |
| move | codes\realtest_code\munich_rent\codes\save_munich_data.R | 05_munich_rent\hpc\codes\save_munich_data.R | BCCD86E44D3BBD2E29EA832999988602 | BCCD86E44D3BBD2E29EA832999988602 |
| move | codes\realtest_code\munich_rent\submit_munich_rent.sub | 05_munich_rent\hpc\submit_munich_rent.sub | 4C4672E9BE1329F86A8AD7E2A46105FD | 4C4672E9BE1329F86A8AD7E2A46105FD |
| move | codes\realtest_code\munich_rent\submit_munich_rent_n1000.sub | 05_munich_rent\hpc\submit_munich_rent_n1000.sub | 4E2C64CF6C2DDFF5B41347CFB4E5F65B | 4E2C64CF6C2DDFF5B41347CFB4E5F65B |
| move | codes\run_codes\munich\baby_bayes.R | 05_munich_rent\local\baby_bayes.R | 8EDC5C1BB9E78BAAEA741A7EC2D38C4F | 8EDC5C1BB9E78BAAEA741A7EC2D38C4F |
| move | codes\run_codes\munich\fit_spatial_reml.R | 05_munich_rent\local\fit_spatial_reml.R | 0F803D62A7CE9448CF6591977507AD78 | 0F803D62A7CE9448CF6591977507AD78 |
| move | codes\run_codes\munich\gibbs_bayes.R | 05_munich_rent\local\gibbs_bayes.R | B8D2CD8B3532518C6539F6ABD8EFE841 | B8D2CD8B3532518C6539F6ABD8EFE841 |
| move | codes\run_codes\munich\gibbs_interaction.R | 05_munich_rent\local\gibbs_interaction.R | 9C0A1288E0E8899AE1ED15FC4BF63AFC | 9C0A1288E0E8899AE1ED15FC4BF63AFC |
| move | codes\run_codes\munich\gibbs_stage_c_full.R | 05_munich_rent\local\gibbs_stage_c_full.R | D22348BA7CD55F6DC3459BC08E7F387C | D22348BA7CD55F6DC3459BC08E7F387C |
| move | codes\run_codes\munich\ls_basis.R | 05_munich_rent\local\ls_basis.R | 6A3132B2720451E7C7B53029CE5C0EFE | 6A3132B2720451E7C7B53029CE5C0EFE |
| move | codes\run_codes\munich\ls_interaction.R | 05_munich_rent\local\ls_interaction.R | 4AE41D9C80024EF4A3732409BE917D47 | 4AE41D9C80024EF4A3732409BE917D47 |
| move | codes\run_codes\munich\ls_interaction_core.cpp | 05_munich_rent\local\ls_interaction_core.cpp | 81F931C45DE492857D5DE0CE4FEB00E3 | 81F931C45DE492857D5DE0CE4FEB00E3 |
| move | codes\run_codes\munich\marginal_utils.R | 05_munich_rent\local\marginal_utils.R | 04498862A8B2131139691A94630EB827 | 04498862A8B2131139691A94630EB827 |
| move | codes\run_codes\munich\plot_modelA_checkpoint.R | 05_munich_rent\local\plot_modelA_checkpoint.R | 19CC7B81FC9E3A2A35A4F6225D9B941C | 19CC7B81FC9E3A2A35A4F6225D9B941C |
| move | codes\run_codes\munich\run_munich_rent_local.R | 05_munich_rent\local\run_munich_rent_local.R | F6441ADF2B6652F379E0A20F0A3EDC3F | F6441ADF2B6652F379E0A20F0A3EDC3F |
| move | codes\run_codes\munich\spatial_utils.R | 05_munich_rent\local\spatial_utils.R | F0216DB4FD9FD76269EDC0F3E897CD40 | F0216DB4FD9FD76269EDC0F3E897CD40 |
| move-bulk-verify | codes\realtest_code\munich_rent | 05_munich_rent\hpc | count=24;bytes=442310 | count=24;bytes=442310 |
| move-bulk-verify | codes\run_codes\munich | 05_munich_rent\local | count=34;bytes=95827427 | count=34;bytes=95827427 |

## docs

| action | source | destination | MD5-before | MD5-after |
|---|---|---|---|---|
| move | baby_bayes_guide.md | docs\baby_bayes_guide.md | AE4082983D41C3938E46C9D67D1A1B46 | AE4082983D41C3938E46C9D67D1A1B46 |
| move | Bayesian_Model_and_Simulation_Sections.docx | docs\Bayesian_Model_and_Simulation_Sections.docx | 2F5808FF8205CA5151B45AC7506014CA | 2F5808FF8205CA5151B45AC7506014CA |
| move | interaction\interaction constract\ls_interaction_bayes_and_experiment.md | docs\interaction_notes\ls_interaction_bayes_and_experiment.md | 5D152F4BC80F14A19D9148650BAABC0C | 5D152F4BC80F14A19D9148650BAABC0C |
| move | interaction\interaction constract\ls_interaction_construction.md | docs\interaction_notes\ls_interaction_construction.md | A3ED7F0A3FA2CAF35EA917CC45B918B7 | A3ED7F0A3FA2CAF35EA917CC45B918B7 |
| move | presenting files\idad_simple_prep.md | docs\presenting\idad_simple_prep.md | 700389EA19B959702314461C495592BF | 700389EA19B959702314461C495592BF |
| move | LS-Spline-Dissertation\00 Maps of Content\MOC - Dissertation.md | docs\vault\00 Maps of Content\MOC - Dissertation.md | 9A0780874B30C6B376CC4F0377BEE529 | 9A0780874B30C6B376CC4F0377BEE529 |
| move | LS-Spline-Dissertation\00 Maps of Content\MOC - Literature.md | docs\vault\00 Maps of Content\MOC - Literature.md | EF09D531665FD93BE9F1A579B5C27DAE | EF09D531665FD93BE9F1A579B5C27DAE |
| move | LS-Spline-Dissertation\00 Maps of Content\MOC - Methods.md | docs\vault\00 Maps of Content\MOC - Methods.md | 1737092A659BD74EB789CBA8C3E28B52 | 1737092A659BD74EB789CBA8C3E28B52 |
| move | LS-Spline-Dissertation\00 Maps of Content\MOC - Simulations.md | docs\vault\00 Maps of Content\MOC - Simulations.md | 73F4AD65E565D2A4FB62D757D6A5D596 | 73F4AD65E565D2A4FB62D757D6A5D596 |
| move | LS-Spline-Dissertation\00 Maps of Content\STUB - Notes still to write.md | docs\vault\00 Maps of Content\STUB - Notes still to write.md | 321E3696F0C299286C4FD56FE924C2B0 | 321E3696F0C299286C4FD56FE924C2B0 |
| move | LS-Spline-Dissertation\10 Methods\ANOVA identifiability.md | docs\vault\10 Methods\ANOVA identifiability.md | 0E3DDF9B646B33068F338FAB0E1B8DA2 | 0E3DDF9B646B33068F338FAB0E1B8DA2 |
| move | LS-Spline-Dissertation\10 Methods\Collapsed Gibbs sampler.md | docs\vault\10 Methods\Collapsed Gibbs sampler.md | 1FB4DAC3ABB711D5619952D94C537072 | 1FB4DAC3ABB711D5619952D94C537072 |
| move | LS-Spline-Dissertation\10 Methods\Kronecker-sum RW2 penalty.md | docs\vault\10 Methods\Kronecker-sum RW2 penalty.md | 013FF4692ABC196119FE15457A3B3F09 | 013FF4692ABC196119FE15457A3B3F09 |
| move | LS-Spline-Dissertation\10 Methods\LS basis.md | docs\vault\10 Methods\LS basis.md | 0370FF9A542BDA098545C0D91F95CB94 | 0370FF9A542BDA098545C0D91F95CB94 |
| move | LS-Spline-Dissertation\10 Methods\Mat├⌐rn covariance function.md | docs\vault\10 Methods\Mat├⌐rn covariance function.md | 4B70D53BC759081D1F4865A6BAF4AE6E | 4B70D53BC759081D1F4865A6BAF4AE6E |
| move | LS-Spline-Dissertation\10 Methods\Mat├⌐rn GP random effect.md | docs\vault\10 Methods\Mat├⌐rn GP random effect.md | 84D724DF0516EC1D8FB8ED1B5C172469 | 84D724DF0516EC1D8FB8ED1B5C172469 |
| move | LS-Spline-Dissertation\10 Methods\RW2 prior.md | docs\vault\10 Methods\RW2 prior.md | AAF68E24D7DA691277842E3B7DD27F0A | AAF68E24D7DA691277842E3B7DD27F0A |
| move | LS-Spline-Dissertation\10 Methods\Tensor-product LS basis.md | docs\vault\10 Methods\Tensor-product LS basis.md | BDBEB0C3FFDC6BB9835C8F86A1DBF260 | BDBEB0C3FFDC6BB9835C8F86A1DBF260 |
| move | LS-Spline-Dissertation\20 Simulations\Interaction 2x2 v2.md | docs\vault\20 Simulations\Interaction 2x2 v2.md | F0C04AF8FD4B332813E26168EAF310CF | F0C04AF8FD4B332813E26168EAF310CF |
| move | LS-Spline-Dissertation\20 Simulations\Prior sensitivity B1-B5.md | docs\vault\20 Simulations\Prior sensitivity B1-B5.md | C428F6C0BBFF05FF9F8B57D2308160D9 | C428F6C0BBFF05FF9F8B57D2308160D9 |
| move | LS-Spline-Dissertation\20 Simulations\REML vs Bayes.md | docs\vault\20 Simulations\REML vs Bayes.md | F4761A9FC5B34E362C6CE97A07E5EA10 | F4761A9FC5B34E362C6CE97A07E5EA10 |
| move | LS-Spline-Dissertation\20 Simulations\Sample size comparison.md | docs\vault\20 Simulations\Sample size comparison.md | B6C89E4D7FFF20B564C6F622023AE479 | B6C89E4D7FFF20B564C6F622023AE479 |
| move | LS-Spline-Dissertation\20 Simulations\Sim1 n1000.md | docs\vault\20 Simulations\Sim1 n1000.md | A3E6415B2AE6B4DC47181FAB518F5588 | A3E6415B2AE6B4DC47181FAB518F5588 |
| move | LS-Spline-Dissertation\30 Applications\Munich rent data.md | docs\vault\30 Applications\Munich rent data.md | 5502C2A016EB0C41FDAF8EB7FBE854A2 | 5502C2A016EB0C41FDAF8EB7FBE854A2 |
| move | LS-Spline-Dissertation\30 Applications\SEER lung cancer queued.md | docs\vault\30 Applications\SEER lung cancer queued.md | 4E6EE11E1F1BB1AF1520C29534206527 | 4E6EE11E1F1BB1AF1520C29534206527 |
| move | LS-Spline-Dissertation\40 Literature\Chib and Greenberg 2010.md | docs\vault\40 Literature\Chib and Greenberg 2010.md | F67A4A1F15A83CCC344C304446E9EFDB | F67A4A1F15A83CCC344C304446E9EFDB |
| move | LS-Spline-Dissertation\40 Literature\Eilers and Marx 2003.md | docs\vault\40 Literature\Eilers and Marx 2003.md | 31D42F50DF56BE78C3F2B9B97C6B3708 | 31D42F50DF56BE78C3F2B9B97C6B3708 |
| move | LS-Spline-Dissertation\40 Literature\Lang and Brezger 2004.md | docs\vault\40 Literature\Lang and Brezger 2004.md | 0B17DD4DD7027C3983B92249F6526CC5 | 0B17DD4DD7027C3983B92249F6526CC5 |
| move | LS-Spline-Dissertation\40 Literature\Nandy Lim Maiti 2017.md | docs\vault\40 Literature\Nandy Lim Maiti 2017.md | 1DE12FCD552012F414B77DDBC10D7CAF | 1DE12FCD552012F414B77DDBC10D7CAF |
| move | LS-Spline-Dissertation\40 Literature\Vehtari et al 2017.md | docs\vault\40 Literature\Vehtari et al 2017.md | CAB0E825C991C97A5C6A7327575DAEFE | CAB0E825C991C97A5C6A7327575DAEFE |
| move | LS-Spline-Dissertation\50 Decisions and Insights\Centering for ANOVA identifiability.md | docs\vault\50 Decisions and Insights\Centering for ANOVA identifiability.md | 3A150481A6FEAE23185F0859D4F96895 | 3A150481A6FEAE23185F0859D4F96895 |
| move | LS-Spline-Dissertation\50 Decisions and Insights\INLA as computational benchmark.md | docs\vault\50 Decisions and Insights\INLA as computational benchmark.md | C9981D3FDF33F8DCDAE6DAA379CE0BCE | C9981D3FDF33F8DCDAE6DAA379CE0BCE |
| move | LS-Spline-Dissertation\50 Decisions and Insights\Novelty stack.md | docs\vault\50 Decisions and Insights\Novelty stack.md | EC382785C9CD8709157E342FAAF072F6 | EC382785C9CD8709157E342FAAF072F6 |
| move | LS-Spline-Dissertation\50 Decisions and Insights\Prior sensitivity is a strength.md | docs\vault\50 Decisions and Insights\Prior sensitivity is a strength.md | A4C0097CFD526D837AD4AAE45B2E2D93 | A4C0097CFD526D837AD4AAE45B2E2D93 |
| move | LS-Spline-Dissertation\50 Decisions and Insights\Why collapsed Gibbs.md | docs\vault\50 Decisions and Insights\Why collapsed Gibbs.md | 6F3CF7245DDB53E038FA092DF9C9FC69 | 6F3CF7245DDB53E038FA092DF9C9FC69 |
| move | LS-Spline-Dissertation\50 Decisions and Insights\Why LS over P-splines.md | docs\vault\50 Decisions and Insights\Why LS over P-splines.md | 8B575341846DF3CC0B002E94115180C1 | 8B575341846DF3CC0B002E94115180C1 |
| move | LS-Spline-Dissertation\50 Decisions and Insights\X1 undercoverage is basis resolution.md | docs\vault\50 Decisions and Insights\X1 undercoverage is basis resolution.md | 3DFF0158BAC3476A0FC4C793E80BEB12 | 3DFF0158BAC3476A0FC4C793E80BEB12 |
| move | LS-Spline-Dissertation\60 Code and HPC\Bug log.md | docs\vault\60 Code and HPC\Bug log.md | 60363C2166497BF30F98E2309AF9DFD8 | 60363C2166497BF30F98E2309AF9DFD8 |
| move | LS-Spline-Dissertation\60 Code and HPC\Hellbender HPC notes.md | docs\vault\60 Code and HPC\Hellbender HPC notes.md | 863881B245A1DEBE363FAE4A63599106 | 863881B245A1DEBE363FAE4A63599106 |
| move | LS-Spline-Dissertation\60 Code and HPC\R file index.md | docs\vault\60 Code and HPC\R file index.md | 795D603AECE2B9AA876216B58F293713 | 795D603AECE2B9AA876216B58F293713 |
| move | LS-Spline-Dissertation\70 Daily\2026-04-28.md | docs\vault\70 Daily\2026-04-28.md | A043502F9AA1D7923A57A0903917A8B1 | A043502F9AA1D7923A57A0903917A8B1 |
| move | LS-Spline-Dissertation\80 Templates\Daily Template.md | docs\vault\80 Templates\Daily Template.md | 01C635287F517E10FA5C4F94F8E02849 | 01C635287F517E10FA5C4F94F8E02849 |
| move | LS-Spline-Dissertation\80 Templates\Decision Template.md | docs\vault\80 Templates\Decision Template.md | AECB24F442BD82136B3A4F11EA421D94 | AECB24F442BD82136B3A4F11EA421D94 |
| move | LS-Spline-Dissertation\80 Templates\Literature Template.md | docs\vault\80 Templates\Literature Template.md | 6D43B8D3B3DD2A0FD3621AFC95434C76 | 6D43B8D3B3DD2A0FD3621AFC95434C76 |
| move | LS-Spline-Dissertation\80 Templates\Method Template.md | docs\vault\80 Templates\Method Template.md | D32B2B1FB4215C140424B9E094641744 | D32B2B1FB4215C140424B9E094641744 |
| move | LS-Spline-Dissertation\80 Templates\Simulation Template.md | docs\vault\80 Templates\Simulation Template.md | 29831B5D20B6578906F4990A82F63363 | 29831B5D20B6578906F4990A82F63363 |
| move | LS-Spline-Dissertation\README.md | docs\vault\README.md | F839027444F2AA0681CE1A63DFAE433B | F839027444F2AA0681CE1A63DFAE433B |
| move-bulk-verify | interaction\interaction constract | docs\interaction_notes | count=7;bytes=580924 | count=7;bytes=580924 |
| move-bulk-verify | interaction\paper | docs\literature\interaction_paper | count=7;bytes=10866018 | count=7;bytes=10866018 |
| move-bulk-verify | paper reading | docs\literature\paper_reading | count=14;bytes=55116650 | count=14;bytes=55116650 |
| move-bulk-verify | presenting files | docs\presenting | count=3;bytes=386809 | count=3;bytes=386809 |
| move-bulk-verify | LS-Spline-Dissertation | docs\vault | count=47;bytes=81117 | count=47;bytes=81117 |

## results

| action | source | destination | MD5-before | MD5-after |
|---|---|---|---|---|
| move | Hellbender-12548977.out | results\_job_logs\Hellbender-12548977.out | 503F54E13763B16D698675668FC3513A | 503F54E13763B16D698675668FC3513A |
| move | Hellbender-12549745.out | results\_job_logs\Hellbender-12549745.out | 790908209FA6994DBBE44856B7FF496E | 790908209FA6994DBBE44856B7FF496E |
| move-bulk-verify | 0311_output | results\0311_output | count=84;bytes=1931985 | count=84;bytes=1931985 |
| move-bulk-verify | all_bayes_results_n400 | results\all_bayes_results_n400 | count=17;bytes=276984 | count=17;bytes=276984 |
| move-bulk-verify | bayes_out_gibss | results\bayes_out_gibss | count=22;bytes=255512 | count=22;bytes=255512 |
| move-bulk-verify | reml_resulyts | results\reml_resulyts | count=57;bytes=410397 | count=57;bytes=410397 |

## _reorg

| action | source | destination | MD5-before | MD5-after |
|---|---|---|---|---|
| move | CANONICAL_CHECK.md | _reorg\CANONICAL_CHECK.md | 6278F40B6EB0D18A6EB970410BE7070F | 6278F40B6EB0D18A6EB970410BE7070F |
| move | CODE_AUDIT.md | _reorg\CODE_AUDIT.md | 242C80642A0120F93478D35A835EEE43 | 242C80642A0120F93478D35A835EEE43 |

## archive

| action | source | destination | MD5-before | MD5-after |
|---|---|---|---|---|
| move | codes\interaction_old\comparison_stage1\canonical_schema.R | archive\interaction_old\comparison_stage1\canonical_schema.R | 76EAB373FE11104E7B7AA7A507D5A104 | 76EAB373FE11104E7B7AA7A507D5A104 |
| move | codes\interaction_old\comparison_stage1\diagnose_centering_bug.R | archive\interaction_old\comparison_stage1\diagnose_centering_bug.R | C71D1AC1EE6975879530CDF5675246CF | C71D1AC1EE6975879530CDF5675246CF |
| move | codes\interaction_old\comparison_stage1\diagnose_f12_absorption.R | archive\interaction_old\comparison_stage1\diagnose_f12_absorption.R | 5562A118CCAC855FCC276592BB753898 | 5562A118CCAC855FCC276592BB753898 |
| move | codes\interaction_old\comparison_stage1\fit_mgcv_interaction_wrapper.R | archive\interaction_old\comparison_stage1\fit_mgcv_interaction_wrapper.R | 36AEDB0A1198C5A660755552719E0D5B | 36AEDB0A1198C5A660755552719E0D5B |
| move | codes\interaction_old\comparison_stage1\fit_ours_interaction_wrapper.R | archive\interaction_old\comparison_stage1\fit_ours_interaction_wrapper.R | A3708FCB0D52DEA355D8429B2896500C | A3708FCB0D52DEA355D8429B2896500C |
| move | codes\interaction_old\comparison_stage1\fit_spatial_reml.R | archive\interaction_old\comparison_stage1\fit_spatial_reml.R | 0F803D62A7CE9448CF6591977507AD78 | 0F803D62A7CE9448CF6591977507AD78 |
| move | codes\interaction_old\comparison_stage1\gibbs_interaction.R | archive\interaction_old\comparison_stage1\gibbs_interaction.R | DB8234D7174E0E839FC69F99530814EA | DB8234D7174E0E839FC69F99530814EA |
| move | codes\interaction_old\comparison_stage1\gibbs_stage_c_full.R | archive\interaction_old\comparison_stage1\gibbs_stage_c_full.R | D22348BA7CD55F6DC3459BC08E7F387C | D22348BA7CD55F6DC3459BC08E7F387C |
| move | codes\interaction_old\comparison_stage1\ls_basis.R | archive\interaction_old\comparison_stage1\ls_basis.R | 6A3132B2720451E7C7B53029CE5C0EFE | 6A3132B2720451E7C7B53029CE5C0EFE |
| move | codes\interaction_old\comparison_stage1\ls_interaction.R | archive\interaction_old\comparison_stage1\ls_interaction.R | 5C2E82A46E7CCCFDAF53E803E9A25C04 | 5C2E82A46E7CCCFDAF53E803E9A25C04 |
| move | codes\interaction_old\comparison_stage1\ls_interaction_core.cpp | archive\interaction_old\comparison_stage1\ls_interaction_core.cpp | 81F931C45DE492857D5DE0CE4FEB00E3 | 81F931C45DE492857D5DE0CE4FEB00E3 |
| move | codes\interaction_old\comparison_stage1\marginal_utils.R | archive\interaction_old\comparison_stage1\marginal_utils.R | 04498862A8B2131139691A94630EB827 | 04498862A8B2131139691A94630EB827 |
| move | codes\interaction_old\comparison_stage1\scenario_A\canonical_schema.R | archive\interaction_old\comparison_stage1\scenario_A\canonical_schema.R | 76EAB373FE11104E7B7AA7A507D5A104 | 76EAB373FE11104E7B7AA7A507D5A104 |
| move | codes\interaction_old\comparison_stage1\scenario_A\compute_recovery_metrics.R | archive\interaction_old\comparison_stage1\scenario_A\compute_recovery_metrics.R | 859E82C58A8BFC034922F344A7B4F0D6 | 859E82C58A8BFC034922F344A7B4F0D6 |
| move | codes\interaction_old\comparison_stage1\scenario_A\fit_inla_wrapper.R | archive\interaction_old\comparison_stage1\scenario_A\fit_inla_wrapper.R | 2B0B8A4DBE48A2999BEDA7C8B20D7FD5 | 2B0B8A4DBE48A2999BEDA7C8B20D7FD5 |
| move | codes\interaction_old\comparison_stage1\scenario_A\fit_mgcv_wrapper.R | archive\interaction_old\comparison_stage1\scenario_A\fit_mgcv_wrapper.R | C05F46CA60878A2F272255F5B13D9787 | C05F46CA60878A2F272255F5B13D9787 |
| move | codes\interaction_old\comparison_stage1\scenario_A\fit_ours_wrapper.R | archive\interaction_old\comparison_stage1\scenario_A\fit_ours_wrapper.R | DB67A88AA37A481DD696F04905F293A7 | DB67A88AA37A481DD696F04905F293A7 |
| move | codes\interaction_old\comparison_stage1\scenario_A\fit_spatial_reml.R | archive\interaction_old\comparison_stage1\scenario_A\fit_spatial_reml.R | 0F803D62A7CE9448CF6591977507AD78 | 0F803D62A7CE9448CF6591977507AD78 |
| move | codes\interaction_old\comparison_stage1\scenario_A\gibbs_stage_c_full.R | archive\interaction_old\comparison_stage1\scenario_A\gibbs_stage_c_full.R | D22348BA7CD55F6DC3459BC08E7F387C | D22348BA7CD55F6DC3459BC08E7F387C |
| move | codes\interaction_old\comparison_stage1\scenario_A\ls_basis.R | archive\interaction_old\comparison_stage1\scenario_A\ls_basis.R | 6A3132B2720451E7C7B53029CE5C0EFE | 6A3132B2720451E7C7B53029CE5C0EFE |
| move | codes\interaction_old\comparison_stage1\scenario_A\marginal_utils.R | archive\interaction_old\comparison_stage1\scenario_A\marginal_utils.R | 04498862A8B2131139691A94630EB827 | 04498862A8B2131139691A94630EB827 |
| move | codes\interaction_old\comparison_stage1\scenario_A\run_stage1_scenarioA_multiseed.R | archive\interaction_old\comparison_stage1\scenario_A\run_stage1_scenarioA_multiseed.R | 51E212C879591ECF140EBB1A22A5381C | 51E212C879591ECF140EBB1A22A5381C |
| move | codes\interaction_old\comparison_stage1\scenario_A\scenarioA_settings.md | archive\interaction_old\comparison_stage1\scenario_A\scenarioA_settings.md | F9F885B826B1FAE8581FBAC3446FB19C | F9F885B826B1FAE8581FBAC3446FB19C |
| move | codes\interaction_old\comparison_stage1\scenario_A\simulate_scenario_A.R | archive\interaction_old\comparison_stage1\scenario_A\simulate_scenario_A.R | BB45788287D6CE2272AEFA40D3706E80 | BB45788287D6CE2272AEFA40D3706E80 |
| move | codes\interaction_old\comparison_stage1\scenario_A\spatial_utils.R | archive\interaction_old\comparison_stage1\scenario_A\spatial_utils.R | F0216DB4FD9FD76269EDC0F3E897CD40 | F0216DB4FD9FD76269EDC0F3E897CD40 |
| move | codes\interaction_old\comparison_stage1\scenario_A\test_fit_inla.R | archive\interaction_old\comparison_stage1\scenario_A\test_fit_inla.R | 2A027C9C67531F87ECD86384768D3426 | 2A027C9C67531F87ECD86384768D3426 |
| move | codes\interaction_old\comparison_stage1\scenario_A\test_fit_mgcv.R | archive\interaction_old\comparison_stage1\scenario_A\test_fit_mgcv.R | 53D3FF54D8F5C0FDF487935FA1C4119D | 53D3FF54D8F5C0FDF487935FA1C4119D |
| move | codes\interaction_old\comparison_stage1\scenario_A\test_fit_ours.R | archive\interaction_old\comparison_stage1\scenario_A\test_fit_ours.R | 833971C6EBC87DBA1A3754832A6103D7 | 833971C6EBC87DBA1A3754832A6103D7 |
| move | codes\interaction_old\comparison_stage1\scenario_A\test_simulate_scenario_A.R | archive\interaction_old\comparison_stage1\scenario_A\test_simulate_scenario_A.R | 76731691C0A6380C4F42547C7D6B4434 | 76731691C0A6380C4F42547C7D6B4434 |
| move | codes\interaction_old\comparison_stage1\simulate_scenario_A.R | archive\interaction_old\comparison_stage1\simulate_scenario_A.R | BB45788287D6CE2272AEFA40D3706E80 | BB45788287D6CE2272AEFA40D3706E80 |
| move | codes\interaction_old\comparison_stage1\simulate_scenario_B.R | archive\interaction_old\comparison_stage1\simulate_scenario_B.R | 45BE1ACA0C8431684B8E442726D582CC | 45BE1ACA0C8431684B8E442726D582CC |
| move | codes\interaction_old\comparison_stage1\spatial_utils.R | archive\interaction_old\comparison_stage1\spatial_utils.R | F0216DB4FD9FD76269EDC0F3E897CD40 | F0216DB4FD9FD76269EDC0F3E897CD40 |
| move | codes\interaction_old\comparison_stage1\test_b_smooth_loose.R | archive\interaction_old\comparison_stage1\test_b_smooth_loose.R | F63E10879E938A2BBF19CFF48B94C590 | F63E10879E938A2BBF19CFF48B94C590 |
| move | codes\interaction_old\comparison_stage1\test_fit_mgcv_interaction.R | archive\interaction_old\comparison_stage1\test_fit_mgcv_interaction.R | ECAC3D845ACB4CB85FEF5044C693D9DC | ECAC3D845ACB4CB85FEF5044C693D9DC |
| move | codes\interaction_old\comparison_stage1\test_fit_ours_interaction.R | archive\interaction_old\comparison_stage1\test_fit_ours_interaction.R | 56A2774E145D3D39A73DD73E5928BFC6 | 56A2774E145D3D39A73DD73E5928BFC6 |
| move | codes\interaction_old\comparison_stage1\test_fit_ours_interaction_ablation.R | archive\interaction_old\comparison_stage1\test_fit_ours_interaction_ablation.R | 09967F6A2FC66D81BB47D457CC5ADAFE | 09967F6A2FC66D81BB47D457CC5ADAFE |
| move | codes\interaction_old\comparison_stage1\test_interaction_on_additive_data.R | archive\interaction_old\comparison_stage1\test_interaction_on_additive_data.R | 865CA3797CFB83732D7AD96F89EFE38F | 865CA3797CFB83732D7AD96F89EFE38F |
| move | codes\interaction_old\comparison_stage1\test_simulate_scenario_B.R | archive\interaction_old\comparison_stage1\test_simulate_scenario_B.R | 793C100CBB284E05430A602FC25E98A8 | 793C100CBB284E05430A602FC25E98A8 |
| move | codes\interaction_old\interaction_week1\gibbs_interaction.R | archive\interaction_old\interaction_week1\gibbs_interaction.R | DB8234D7174E0E839FC69F99530814EA | DB8234D7174E0E839FC69F99530814EA |
| move | codes\interaction_old\interaction_week1\ls_interaction.R | archive\interaction_old\interaction_week1\ls_interaction.R | 5C2E82A46E7CCCFDAF53E803E9A25C04 | 5C2E82A46E7CCCFDAF53E803E9A25C04 |
| move | codes\interaction_old\interaction_week1\ls_interaction_core.cpp | archive\interaction_old\interaction_week1\ls_interaction_core.cpp | 81F931C45DE492857D5DE0CE4FEB00E3 | 81F931C45DE492857D5DE0CE4FEB00E3 |
| move | codes\interaction_old\interaction_week1\run_interaction_2x2.R | archive\interaction_old\interaction_week1\run_interaction_2x2.R | FEC2DFA6B20E4D74BE2B1B1B62D15210 | FEC2DFA6B20E4D74BE2B1B1B62D15210 |
| move | codes\interaction_old\interaction_week1\run_interaction_2x2_p4.R | archive\interaction_old\interaction_week1\run_interaction_2x2_p4.R | AD7482239FB7C92D40CE804B52695F8F | AD7482239FB7C92D40CE804B52695F8F |
| move | codes\interaction_old\interaction_week1\run_interaction_2x2_v2.R | archive\interaction_old\interaction_week1\run_interaction_2x2_v2.R | E471F378EDC94FC2CD3E17DA7F5C49D6 | E471F378EDC94FC2CD3E17DA7F5C49D6 |
| move | codes\interaction_old\interaction_week1\run_interaction_null_test.R | archive\interaction_old\interaction_week1\run_interaction_null_test.R | 7F9CEB01483AAAADDD1638B9A938AD14 | 7F9CEB01483AAAADDD1638B9A938AD14 |
| move | codes\interaction_old\interaction_week1\scenarioA_settings.md | archive\interaction_old\interaction_week1\scenarioA_settings.md | 30DB4F75AA14243D63701316C58602C4 | 30DB4F75AA14243D63701316C58602C4 |
| move | codes\poisson_extension\gibbs_stage_c_poisson.R | archive\poisson_extension\gibbs_stage_c_poisson.R | 8F5DA9FDCD8CBCF2AB04245D7094C411 | 8F5DA9FDCD8CBCF2AB04245D7094C411 |
| move | codes\poisson_extension\michigan_lung_data\run_michigan_lung_poisson.R | archive\poisson_extension\michigan_lung_data\run_michigan_lung_poisson.R | 3A74CBDB4909FE1C2A3E2B83A6AEB812 | 3A74CBDB4909FE1C2A3E2B83A6AEB812 |
| move | codes\poisson_extension\michigan_lung_data\run_michigan_lung_poisson.sub | archive\poisson_extension\michigan_lung_data\run_michigan_lung_poisson.sub | B2802DA8794BF2526B9A9A578C377B24 | B2802DA8794BF2526B9A9A578C377B24 |
| move | codes\poisson_extension\sim1_2_code\run_sim_poisson_1.R | archive\poisson_extension\sim1_2_code\run_sim_poisson_1.R | 3E2246789D43EE0745389BFB75E18758 | 3E2246789D43EE0745389BFB75E18758 |
| move | codes\poisson_extension\sim1_2_code\run_sim_poisson_1.sub | archive\poisson_extension\sim1_2_code\run_sim_poisson_1.sub | 04F09F53565760382D8DB6585C36E907 | 04F09F53565760382D8DB6585C36E907 |
| move | codes\poisson_extension\sim1_2_code\run_sim_poisson_2.R | archive\poisson_extension\sim1_2_code\run_sim_poisson_2.R | 0D1D19B6B9D4D83239D4539FC4AB2C59 | 0D1D19B6B9D4D83239D4539FC4AB2C59 |
| move | codes\poisson_extension\sim1_2_code\run_sim_poisson_2.sub | archive\poisson_extension\sim1_2_code\run_sim_poisson_2.sub | 979BBE7D1B0FE30FD06F4AE14C8B3692 | 979BBE7D1B0FE30FD06F4AE14C8B3692 |
| move-bulk-verify | codes\interaction_old | archive\interaction_old | count=127;bytes=827790924 | count=127;bytes=827790924 |
| move-bulk-verify | codes\poisson_extension | archive\poisson_extension | count=26;bytes=610663 | count=26;bytes=610663 |

## archive/_to_delete/

| action | source | destination | MD5-before | MD5-after |
|---|---|---|---|---|
| park | 03_bayes_main\n1000\.Rhistory | archive\_to_delete\03_bayes_main__n1000__.Rhistory | D41D8CD98F00B204E9800998ECF8427E | D41D8CD98F00B204E9800998ECF8427E |
| park | 03_bayes_main\n1000\bayes_overlay\.Rhistory | archive\_to_delete\03_bayes_main__n1000__bayes_overlay__.Rhistory | D41D8CD98F00B204E9800998ECF8427E | D41D8CD98F00B204E9800998ECF8427E |
| park | 03_bayes_main\n1000\bayes_test_sim1\.Rhistory | archive\_to_delete\03_bayes_main__n1000__bayes_test_sim1__.Rhistory | D41D8CD98F00B204E9800998ECF8427E | D41D8CD98F00B204E9800998ECF8427E |
| park | 05_munich_rent\hpc\.Rhistory | archive\_to_delete\05_munich_rent__hpc__.Rhistory | D41D8CD98F00B204E9800998ECF8427E | D41D8CD98F00B204E9800998ECF8427E |
| park | codes\.RData | archive\_to_delete\codes__.RData | B727EEEDD8A1A2B05C3C12926758BA3D | B727EEEDD8A1A2B05C3C12926758BA3D |
| park | codes\.Rhistory | archive\_to_delete\codes__.Rhistory | AF9D71EAF4C754F82E314724BA97BA53 | AF9D71EAF4C754F82E314724BA97BA53 |
| park | Cubic_Spline_with_Spatial_Data_april.pdf | archive\_to_delete\Cubic_Spline_with_Spatial_Data_april.pdf | B2FBA11DFEA4063E7E1804ACF0C11F6F | B2FBA11DFEA4063E7E1804ACF0C11F6F |
| park | Cubic_Spline_with_Spatial_Data_Feb25.pdf | archive\_to_delete\Cubic_Spline_with_Spatial_Data_Feb25.pdf | 9198FDAADD665E45ACEAC464773762C5 | 9198FDAADD665E45ACEAC464773762C5 |
| park | .RData | archive\_to_delete\root__.RData | D9189C6E5105DC2A5800503B2CDFEC88 | D9189C6E5105DC2A5800503B2CDFEC88 |
| park | .Rhistory | archive\_to_delete\root__.Rhistory | 0515E4B2572CA4B04015C738DEDEDB5A | 0515E4B2572CA4B04015C738DEDEDB5A |

## Self-check

**(a) Byte preservation â€” every MD5-before == MD5-after:**

Zero mismatches. All file actions preserved bytes exactly (verified inline at action time in stages Aâ€“D, and re-hashed at destination here).

**(b) Propagation â€” every propagated core file's MD5 == its `lib/` canonical:**

Checked every code file named like one of the 12 canonical files, excluding `archive/` (frozen stale snapshot), `06_medsat_pilot/` (untouched), `lib/` (the canonical itself), and `_reorg/`.

Zero mismatches. Every propagated/moved copy of a canonical file (in `02_reml`, `03_bayes_main`, `04_interaction`, `05_munich_rent`) is byte-identical to its `lib/` source.

## Counts

- Files moved: 220 individual code-file move rows (of which 53 into `archive/`), plus 19 whole-directory bulk-verify rows covering all binary/data content.
- Files copied: 61 (copy-to-lib: 12, copy-from-lib: 38, copy-as-is: 11).
- Files parked (removal candidates, NOT deleted): 10.
- Files deleted: 0.

> Note: `archive/poisson_extension/michigan_lung_data/.RData` and `.Rhistory` were deliberately **left frozen** inside the archived snapshot (not parked). Their presence in `archive/` rather than `archive/_to_delete/` is itself the record of that decision.

