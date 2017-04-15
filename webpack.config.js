var webpack = require("webpack");
// var CopyWebpackPlugin = require('copy-webpack-plugin');
module.exports = {
  entry: {
    application:'./index.js',
    worker: './GeneratorWorker.js'
  },
  target: 'node',
  devtool: 'source-map',
  module: {
    loaders: [
      {
        test: /\.js$/,
        loader: 'babel',
        exclude: /node_modules/,
        babelrc:false,
        query: {
          presets: ['es2015', 'stage-0'],
        }

      },
      {test: /\.json$/, loader: "json"},
      {test: /\.ts$/, loader: "ts"},
    ]
  },
  devServer: {
    headers: { "Access-Control-Allow-Origin": "*" }
  },
  output: {
    path: __dirname + '/build/app/',
    publicPath: "/app/",
    filename: '[name].js',
  }
};

