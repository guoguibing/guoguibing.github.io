// Copyright (C) 2014 Guibing Guo
//
// This file is part of LibRec.
//
// LibRec is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// LibRec is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with LibRec. If not, see <http://www.gnu.org/licenses/>.
//

package librec.ranking;

import java.util.ArrayList;
import java.util.List;

import librec.data.Configuration;
import librec.data.DenseMatrix;
import librec.data.DenseVector;
import librec.data.SparseMatrix;
import librec.intf.SocialRecommender;
import librec.util.Randoms;
import librec.util.Strings;

/**
 * FST: Factored Similarity Model with Trust for Item Recommendation
 * 
 * Under development...
 * 
 * @author guoguibing
 * 
 */
@Configuration("binThold, rho, alpha, z, factors, lRate, maxLRate, regU, regI, regB, iters")
public class FST extends SocialRecommender {

	private int rho;
	private float alpha, beta, z;

	private DenseMatrix X, Y;

	public FST(SparseMatrix trainMatrix, SparseMatrix testMatrix, int fold) {
		super(trainMatrix, testMatrix, fold);

		isRankingPred = true;
	}

	@Override
	protected void initModel() throws Exception {
		P = new DenseMatrix(numUsers, numFactors);
		Q = new DenseMatrix(numUsers, numFactors);

		X = new DenseMatrix(numItems, numFactors);
		Y = new DenseMatrix(numItems, numFactors);

		P.init(smallValue);
		Q.init(smallValue);

		X.init(smallValue);
		Y.init(smallValue);

		itemBias = new DenseVector(numItems);
		itemBias.init(smallValue);

		rho = algoOptions.getInt("-rho");
		alpha = algoOptions.getFloat("-alpha");
		beta = algoOptions.getFloat("-beta");
		z = algoOptions.getFloat("-z");

		userItemsCache = trainMatrix.rowColumnsCache(cacheSpec);
		itemUsersCache = trainMatrix.columnRowsCache(cacheSpec);
		userFriendsCache = socialMatrix.rowColumnsCache(cacheSpec);
	}

	@Override
	protected void buildModel() throws Exception {

		for (int iter = 1; iter <= numIters; iter++) {

			loss = 0;

			DenseMatrix PS = new DenseMatrix(numUsers, numFactors);
			DenseMatrix QS = new DenseMatrix(numUsers, numFactors);

			DenseMatrix XS = new DenseMatrix(numItems, numFactors);
			DenseMatrix YS = new DenseMatrix(numItems, numFactors);

			// update throughout each user-item-rating (u, j, ruj) cell
			for (int u : trainMatrix.rows()) {
				List<Integer> ratedItems = userItemsCache.get(u);
				List<Integer> trustUsers = userFriendsCache.get(u);
				double wt = trustUsers.size() > 0 ? Math.pow(trustUsers.size(), -z) : 0;

				for (int i : ratedItems) {
					double rui = trainMatrix.get(u, i);

					// sample a set of items unrated by user u
					List<Integer> js = new ArrayList<>();
					int len = 0;
					while (len < rho) {
						int j = Randoms.uniform(numItems);
						if (ratedItems.contains(j) || js.contains(j))
							continue;

						js.add(j);
						len++;
					}

					// user similarity
					double sum_vi = 0;
					double[] sum_vif = new double[numFactors];
					int cnt_v = 0;
					List<Integer> ratingUsers = itemUsersCache.get(i);
					for (int v : ratingUsers) {
						if (u != v) {
							sum_vi += DenseMatrix.rowMult(P, v, Q, u);
							cnt_v++;

							for (int f = 0; f < numFactors; f++) {
								sum_vif[f] += P.get(v, f);
							}
						}
					}
					double w_vi = cnt_v > 0 ? Math.pow(cnt_v, -beta) : 0;

					// item similarity
					double sum_wi = 0;
					int cnt_w = 0;
					double[] sum_wif = new double[numFactors];
					for (int w : ratedItems) {
						if (w != i) {
							sum_wi += DenseMatrix.rowMult(X, w, Y, i);
							cnt_w++;

							for (int f = 0; f < numFactors; f++) {
								sum_wif[f] += X.get(w, f);
							}
						}
					}
					double w_wi = cnt_w > 0 ? Math.pow(cnt_w, -alpha) : 0;
					double w_wj = ratedItems.size() > 0 ? Math.pow(ratedItems.size(), -alpha) : 0;

					// trust influence
					double sum_ti = 0;
					double[] sum_tf = new double[numFactors];
					for (int t : trustUsers) {
						sum_ti += DenseMatrix.rowMult(P, t, Y, i);

						for (int f = 0; f < numFactors; f++) {
							sum_tf[f] += P.get(t, f);
						}
					}

					// update for each item j unrated by user u
					double[] xs = new double[numFactors];
					double[] ws = new double[numFactors];
					double[] ys = new double[numFactors];
					for (int j : js) {

						List<Integer> Cj = itemUsersCache.get(j);
						double sum_vj = 0;
						double[] sum_vjf = new double[numFactors];
						int cnt_j = 0;
						for (int v : Cj) {
							sum_vj += DenseMatrix.rowMult(P, v, Q, u);
							cnt_j++;

							for (int f = 0; f < numFactors; f++) {
								sum_vjf[f] += P.get(v, f);
							}
						}
						double w_vj = cnt_j > 0 ? Math.pow(cnt_j, -beta) : 0;

						double sum_wj = 0;
						double[] sum_wjf = new double[numFactors];
						for (int w : ratedItems) {
							sum_wj += DenseMatrix.rowMult(X, w, Y, j);

							for (int f = 0; f < numFactors; f++) {
								sum_wjf[f] += X.get(w, f);
							}
						}

						double sum_tj = 0;
						for (int t : trustUsers) {
							sum_tj += DenseMatrix.rowMult(P, t, Y, j);
						}

						double bi = itemBias.get(i), bj = itemBias.get(j);
						double pui = bi + w_vi * sum_vi + w_wi * sum_wi + wt * sum_ti;
						double puj = bj + w_vj * sum_vj + w_wj * sum_wj + wt * sum_tj;
						double ruj = 0;
						double eij = (rui - ruj) - (pui - puj);

						loss += eij * eij;

						// update bi
						itemBias.add(i, -lRate * (-eij + regB * bi));

						// update bj
						itemBias.add(j, -lRate * (eij - regB * bj));

						loss += regB * bi * bi - regB * bj * bj;

						// update quf, yif, yjf
						for (int f = 0; f < numFactors; f++) {
							double quf = Q.get(u, f);
							double yif = Y.get(i, f), yjf = Y.get(j, f);

							double delta = eij * (w_vj * sum_vjf[f] - w_vi * sum_vif[f]) + regU * quf;
							QS.add(u, f, -lRate * delta);

							loss += regU * quf * quf;

							delta = eij * (-w_wi * sum_wif[f] - wt * sum_tf[f]) + regI * yif;
							YS.add(i, f, -lRate * delta);

							delta = eij * (w_wj * sum_wjf[f] + wt * sum_tf[f]) - regI * yjf;
							YS.add(j, f, -lRate * delta);

							loss += regI * yif * yif - regI * yjf * yjf;

							xs[f] += eij * (-w_vi) * quf;
							ws[f] += eij * (w_wj * yjf - w_wi * yif);
							ys[f] += eij * wt * (yjf - yif);
						}

						// update pvf for v in cj
						for (int v : Cj) {
							for (int f = 0; f < numFactors; f++) {
								double pvf = P.get(v, f);
								double delta = eij * w_vj * Q.get(u, f) - regU * pvf;
								PS.add(v, f, -lRate * delta);

								loss -= regU * pvf * pvf;
							}
						}

					}

					// update pvf for v in Ci
					for (int v : ratingUsers) {
						if (v != u) {
							for (int f = 0; f < numFactors; f++) {
								double pvf = P.get(v, f);
								double delta = xs[f] / rho + regU * pvf;
								PS.add(v, f, -lRate * delta);

								loss += regU * pvf * pvf;
							}
						}
					}

					// update xwf for w in Ru
					for (int w : ratedItems) {
						if (w != i) {
							for (int f = 0; f < numFactors; f++) {
								double xwf = X.get(w, f);
								double delta = ws[f] / rho + regI * xwf;
								XS.add(w, f, -lRate * delta);

								loss += regI * xwf * xwf;
							}
						}
					}

					// update qtf for t in Tu
					for (int t : trustUsers) {
						for (int f = 0; f < numFactors; f++) {
							double ptf = P.get(t, f);
							double delta = ys[f] / rho + regU * ptf;
							PS.add(t, f, -lRate * delta);

							loss += regU * ptf * ptf;
						}
					}

				}

			}

			P = P.add(PS);
			Q = Q.add(QS);
			X = X.add(XS);
			Y = Y.add(YS);

			loss *= 0.5;

			if (isConverged(iter))
				break;
		}
	}

	@Override
	public double predict(int u, int i) throws Exception {

		// user similarity
		double sum_c = 0;
		int count = 0;
		List<Integer> ratingUsers = itemUsersCache.get(i);
		for (int v : ratingUsers) {
			if (v != u) {
				sum_c += DenseMatrix.rowMult(P, v, Q, u);
				count++;
			}
		}
		double wc = count > 0 ? Math.pow(count, -beta) : 0;

		// item similarity
		double sum_r = 0;
		count = 0;
		List<Integer> ratedItems = userItemsCache.get(u);
		for (int w : ratedItems) {
			if (w != i) {
				sum_r += DenseMatrix.rowMult(X, w, Y, i);
				count++;
			}
		}
		double wr = count > 0 ? Math.pow(count, -alpha) : 0;

		// trust influence
		double sum_t = 0;
		List<Integer> trustUsers = userFriendsCache.get(u);
		for (int t : trustUsers) {
			sum_t += DenseMatrix.rowMult(P, t, Y, i);
		}
		double wt = trustUsers.size() > 0 ? Math.pow(trustUsers.size(), -z) : 0;

		return itemBias.get(i) + wc * sum_c + wr * sum_r + wt * sum_t;
	}

	@Override
	public String toString() {
		return Strings.toString(new Object[] { binThold, rho, alpha, beta, z, numFactors, initLRate, maxLRate, regU,
				regI, regB, numIters });
	}
}