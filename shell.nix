{ pkgs ? import <nixpkgs> {} }:

pkgs.mkShell {
  buildInputs = with pkgs; [
    entr
    zig
    python312
  ];

  shellHook = ''
    echo "Welcome to the development environment for the colloidal simulation project!"
    echo "Zig version: $(zig version)"
    echo "Python version: $(python --version)"
  '';
}