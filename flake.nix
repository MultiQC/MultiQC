{
  inputs = {
    mach-nix.url = "github:DavHau/mach-nix";
    pypi.url = "github:DavHau/pypi-deps-db";
    pypi.flake = false;
    mach-nix.inputs.pypi-deps-db.follows = "pypi";
  };
  outputs = {self, nixpkgs, ...}@inputs:
    let
      target-system = "x86_64-linux";
      pkgs = import nixpkgs {system = target-system;};

      # this will find all dependencies in setup.py
      mkPackage = system: inputs.mach-nix.lib."${system}".buildPythonPackage {src=./.;
                                                                              python = "python38";
                                                                              requirementsExtra = "setuptools";
                                                                             };
      mkPython = system: inputs.mach-nix.lib."${system}".mkPython {
        python = "python38";
        # extra packages needed for development
        requirements = ''python-lsp-server
                         rich-click
                         black'';

        # hack because black requires tomli (currently broken in mach-nix, issue #484)
        _.black.propagatedBuildInputs.mod = pySelf: self: oldVal: oldVal ++ [ pySelf.tomli ];
        packagesExtra = [ (mkPackage system) ];
      };
    in
  {
    packages.${target-system}.default = mkPackage target-system;
    devShells.${target-system}.default = pkgs.mkShell{
      buildInputs=[(mkPython target-system) pkgs.nodePackages.prettier pkgs.nodePackages.markdownlint-cli];
    };

  };
}
