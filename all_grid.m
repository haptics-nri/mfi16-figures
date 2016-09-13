gs_acc_14v_p = icra17_grid(data14.keys, features14, split_idx14, 20:10:60, {'perceptual'}, .1:.05:.4, .05:.05:.2, 10:20:100);
gs_acc_14v_n = icra17_grid(data14.keys, features14, split_idx14, [3 5 15 25], {'naive'}, 0, .05:.05:.3, 10:20:100);
gs_acc_14b_p = icra17_grid(data14.keys, bfeatures14, bsplit_idx14, 20:10:60, {'perceptual'}, .1:.05:.4, .05:.05:.2, 10:20:100);
gs_acc_14b_n = icra17_grid(data14.keys, bfeatures14, bsplit_idx14, [3 5 15 25], {'naive'}, 0, .05:.05:.3, 10:20:100);
gs_acc_38v_p = icra17_grid(data38.keys, features38, split_idx38, 20:10:60, {'perceptual'}, .1:.05:.4, .05:.05:.2, 10:20:100);
gs_acc_38v_n = icra17_grid(data38.keys, features38, split_idx38, [3 5 15 25], {'naive'}, 0, .05:.05:.3, 10:20:100);
gs_acc_38b_p = icra17_grid(data38.keys, bfeatures38, bsplit_idx38, 20:10:60, {'perceptual'}, .1:.05:.4, .05:.05:.2, 10:20:100);
gs_acc_38b_n = icra17_grid(data38.keys, bfeatures38, bsplit_idx38, [3 5 15 25], {'naive'}, 0, .05:.05:.3, 10:20:100);

%%

icra17_test(data14.keys, features14, split_idx14, gs_acc_14v_p, 20:10:60, {'perceptual'}, .1:.05:.4, .05:.05:.2, 10:20:100);
icra17_test(data14.keys, features14, split_idx14, gs_acc_14v_n, [3 5 15 25], {'naive'}, 0, .05:.05:.3, 10:20:100);
icra17_test(data14.keys, bfeatures14, bsplit_idx14, gs_acc_14b_p, 20:10:60, {'perceptual'}, .1:.05:.4, .05:.05:.2, 10:20:100);
icra17_test(data14.keys, bfeatures14, bsplit_idx14, gs_acc_14b_n, [3 5 15 25], {'naive'}, 0, .05:.05:.3, 10:20:100);
icra17_test(data38.keys, features38, split_idx38, gs_acc_38v_p, 20:10:60, {'perceptual'}, .1:.05:.4, .05:.05:.2, 10:20:100);
icra17_test(data38.keys, features38, split_idx38, gs_acc_38v_n, [3 5 15 25], {'naive'}, 0, .05:.05:.3, 10:20:100);
icra17_test(data38.keys, bfeatures38, bsplit_idx38, gs_acc_38b_p, 20:10:60, {'perceptual'}, .1:.05:.4, .05:.05:.2, 10:20:100);
icra17_test(data38.keys, bfeatures38, bsplit_idx38, gs_acc_38b_n, [3 5 15 25], {'naive'}, 0, .05:.05:.3, 10:20:100);
