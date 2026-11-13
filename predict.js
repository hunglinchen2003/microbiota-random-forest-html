// Prediction formulas copied from the trained randomForest object.
// Split points are little-endian IEEE-754 doubles, base64-encoded,
// so the browser compares the same bits R used.

function b64ToF64(b64) {
  var alpha = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
  var map = {};
  var i;
  for (i = 0; i < alpha.length; i++) map[alpha.charAt(i)] = i;
  var clean = [];
  for (i = 0; i < b64.length; i++) {
    var ch = b64.charAt(i);
    if (ch !== "=" && ch !== "\n" && ch !== "\r") clean.push(ch);
  }
  var n = clean.length;
  var outLen = (n * 3) >> 2;
  var bytes = new Array(outLen);
  var p = 0;
  for (i = 0; i < n; i += 4) {
    var a = map[clean[i]];
    var b = map[clean[i + 1]];
    var c = (i + 2 < n) ? map[clean[i + 2]] : 0;
    var d = (i + 3 < n) ? map[clean[i + 3]] : 0;
    var triple = (a << 18) | (b << 12) | (c << 6) | d;
    if (p < outLen) bytes[p++] = (triple >> 16) & 255;
    if (p < outLen) bytes[p++] = (triple >> 8) & 255;
    if (p < outLen) bytes[p++] = triple & 255;
  }
  var nNum = bytes.length / 8;
  var out = new Array(nNum);
  var buf = new ArrayBuffer(8);
  var view = new DataView(buf);
  var u8 = new Uint8Array(buf);
  var j, k;
  for (j = 0; j < nNum; j++) {
    for (k = 0; k < 8; k++) u8[k] = bytes[j * 8 + k];
    out[j] = view.getFloat64(0, true);
  }
  return out;
}

function hydrateModel(model) {
  var t;
  for (t = 0; t < model.trees.length; t++) {
    if (typeof model.trees[t].S === "string") {
      model.trees[t].S = b64ToF64(model.trees[t].S);
    }
  }
  if (typeof model.importance.raw === "string") {
    model.importance.rawValues = b64ToF64(model.importance.raw);
  }
  return model;
}

// caret::varImp.train with scale = TRUE:
//   shifted = raw - min(raw)
//   scaled  = shifted / max(shifted) * 100
function scaleImportance(raw) {
  var min = raw[0];
  var max = raw[0];
  var i;
  for (i = 1; i < raw.length; i++) {
    if (raw[i] < min) min = raw[i];
    if (raw[i] > max) max = raw[i];
  }
  var shifted = new Array(raw.length);
  var shiftedMax = 0;
  for (i = 0; i < raw.length; i++) {
    shifted[i] = raw[i] - min;
    if (shifted[i] > shiftedMax) shiftedMax = shifted[i];
  }
  var out = new Array(raw.length);
  for (i = 0; i < raw.length; i++) {
    out[i] = shiftedMax === 0 ? 0 : (shifted[i] / shiftedMax) * 100;
  }
  return out;
}

function meanAccuracy(values) {
  var s = 0;
  var i;
  for (i = 0; i < values.length; i++) s += values[i];
  return s / values.length;
}

// classForest rule used by randomForest::predict:
// choose the class that maximises (vote proportion) / cutoff.
// proportion = votes / ntree, and ntree is the same for every class,
// so the comparison is votes[c] / cutoff[c].
// The fitted model uses cutoff = 1/3 for every class, which makes this
// the same as the largest vote count.
// randomForest breaks an exact tie with unif_rand(). This function does
// not draw that random number. It reports tie = true and lists every
// leading class. classIndex is only the deterministic first leader.
function predictCase(model, x) {
  var nclass = model.classes.length;
  var votes = new Array(nclass);
  var c;
  for (c = 0; c < nclass; c++) votes[c] = 0;

  var trees = model.trees;
  var t, node, tree, V, S, L, R, C, f;
  for (t = 0; t < trees.length; t++) {
    tree = trees[t];
    V = tree.V;
    S = tree.S;
    L = tree.L;
    R = tree.R;
    C = tree.C;
    node = 0;
    var guard = 0;
    while (V[node] !== 0) {
      if (++guard > V.length) throw new Error("tree did not reach a leaf");
      f = V[node] - 1;
      if (x[f] <= S[node]) node = L[node];
      else node = R[node];
    }
    votes[C[node] - 1] += 1;
  }

  var cutoff = model.cutoff;
  var equalCutoff = true;
  for (c = 1; c < nclass; c++) {
    if (cutoff[c] !== cutoff[0]) equalCutoff = false;
  }

  var best = 0;
  var leaders = [0];
  if (equalCutoff) {
    for (c = 1; c < nclass; c++) {
      if (votes[c] > votes[best]) {
        best = c;
        leaders = [c];
      } else if (votes[c] === votes[best]) {
        leaders.push(c);
      }
    }
  } else {
    var bestScore = votes[0] / cutoff[0];
    for (c = 1; c < nclass; c++) {
      var score = votes[c] / cutoff[c];
      if (score > bestScore) {
        best = c;
        bestScore = score;
        leaders = [c];
      } else if (score === bestScore) {
        leaders.push(c);
      }
    }
  }

  var ntree = trees.length;
  var prob = new Array(nclass);
  var leaderNames = new Array(leaders.length);
  for (c = 0; c < nclass; c++) prob[c] = votes[c] / ntree;
  for (c = 0; c < leaders.length; c++) leaderNames[c] = model.classes[leaders[c]];

  return {
    classIndex: best,
    species: model.classes[best],
    leaders: leaderNames,
    tie: leaders.length > 1,
    votes: votes,
    prob: prob
  };
}
